#include "lung.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <eigen_3_3_4/Dense>



// Constructor
//**************************************************************/
lung::lung(controlProperties *contProp_, systemProperties *sysProp_){

    pi = atan(1.0)*4.0;

    // Write class members
    //**********************************************************/
    contProp            = contProp_;
    sysProp             = sysProp_;


    // Initialize lung variables and pass values
    //**********************************************************/
    dLimit              = sysProp->dL;
    FRC                 = sysProp->FRC;
    washout             = contProp->washout;
    species             = contProp->species;
    dt                  = contProp->dt;
    Tend                = contProp->Tend;
    nbrBreaths          = contProp->nbrBreaths;


    // Set bifurcation parameters
    maxGenTrumpLob      = sysProp->maxGenLb;
    patm                = sysProp->p0;
    r                   = sysProp->r;
    eta                 = sysProp->eta;
    kappaLb[0]           = pow((1.-r),(1./eta));
    kappaLb[1]           = pow((r),(1./eta));

    // Set homothety ratio in the trumpet lobule (kappa_hat for cross-section)
    kappa               = sysProp->kappa;
    kappa_hat          = 2.*kappa*kappa;

    // Set inlet concentration
    cin                 = contProp->cin;
    bc                  = contProp->bc;

    // Set fluid properties
    mu                  = sysProp->mu;
    rho                 = sysProp->rho;

    // magnification factor for total lobular resistance
    R_fac               = sysProp->R_fac;

    // flag for trumpet lobule scaling
    scalingLbL           = bool(sysProp->scalingLbL);


    // Initialize remaining variables
    //**********************************************************/
    Nt = 1;
    indAbs_Nt = gen_Nt = 0;

    TVm = 0;
    TBm = 0;
    E          = 1e09;

    NntotD = NztotD = NntotLb = NztotLb = 0;

    firstNodeIndD = firstZoneIndD = firstCoordIndD = firstConnIndD = 0;
    firstNodeIndLb = firstZoneIndLb = firstCoordIndLb = firstConnIndLb = 0;

    umin       = 0;

    // Set input flow
    Qin        = 0.;
    Qouttot    = 0.;

    TAWV = TLbV = TLbV_aft_mod = TRLV = TLbVRed = 0;
    totalCompScaling = 0.;
    nbrEndDucts = nbrDucts = nbrEndDuctsCheck = endDuctMinGen = endDuctMinGen = 0;
    nbrLobules = nbrLobulesRed = 0;

    // Initialize Trachea
    //**********************************************************/
    pTrachea = new duct(contProp, sysProp);

    pTrachea->gen = 0;
    pTrachea->ind = 0;
    pTrachea->indAbs = 0;
    pTrachea->endCAW = 0;

    // Scale dimensions with FRC
    pTrachea->d = pow(FRC/0.0035,1./3.)*0.018;
    pTrachea->l = pow(FRC/0.0035,1./3.)*0.122 * 2.; // trachea plus larynx, mounth and apperatus
    pTrachea->lc = pTrachea->l;
    lmin = pTrachea->l;
    dmin = pTrachea->d;

    // For visualization
    pTrachea->x[0] = 0;
    pTrachea->x[1] = pTrachea->l;
    pTrachea->y[0] = 0;
    pTrachea->y[1] = 0;
    pTrachea->z[0] = 0;
    pTrachea->z[1] = 0;
    pTrachea->phi  = 0;
    pTrachea->phi_z= 0;

    // Initialize lobule template
    //**********************************************************/
    pLbTroot = new lobule(contProp, sysProp);

    // Lobule volume in template
    VLbT = 0;

    // Compliance of lobules
    pLbTroot->E = E;

    // For visualization
    Phi[0]  = -pi*1.0/04.;
    Phi[1]  =  pi*1.0/05.;
    Phi[2]  = -pi*1.0/33.;
    Phi[3]  =  pi*1.0/25.;

    // Initialize counter for frames
    frame = 0;
}

// Destructor
//**********************************************************/
lung::~lung(){

    // Delete output data arrays
    delete[] pZoneTypesD;
    delete[] pZoneTypesLb;
    delete[] pConnectivityD;
    delete[] pConnectivityLb;

    delete[] pCoordinatesD;
    delete[] pCoordinatesLb;

    delete[] pZonalVelocityD;
    delete[] pNodalVelocityLb;
    delete[] pNodalPressureD;
    delete[] pNodalConcentrationID;
    delete[] pNodalConcentrationILb;
    delete[] pNodalConcentrationIID;
    delete[] pNodalConcentrationIILb;

    delete[] pDimVarD;
    delete[] pDimVarLb;
    delete[] pCenteringD;
    delete[] pCenteringLb;
    delete[] pVarnamesD;
    delete[] pVarnamesLb;

    delete[] pAllOutputVariablesD;
    delete[] pAllOutputVariablesLb;

    // Create rubber that moves through the duct tree and deletes its branches
    duct* rubberD = new duct(contProp, sysProp);

    // Deconstruct lung
    deconstructLung(pTrachea, rubberD);

    delete rubberD;

}


// Destructor
//**************************************************************/
void lung::deconstructLung(duct* parentDuct, duct* rubberD){

    // point on current duct
    rubberD = parentDuct;

    if (rubberD->reachedEnd(dLimit)){

        /*
        // Create rubber that moves through the lobule tree and delets its branches
        lobule* rubberA = new lobule(contProp, sysProp);

        // Delete trumpet lobule compartement
        rubberA = parentDuct->pLobule;

        delete rubberA;
        */

        return;
    }

    // Call next generation
    deconstructLung(rubberD->pMajorDaughter,rubberD);
    deconstructLung(rubberD->pMinorDaughter,rubberD);

    // Delete current duct
    delete parentDuct;

}

// Generate Morphology of Conducting Airways
//**********************************************************/
void lung::genMorphCondAirways(duct* parentDuct){

    // Set absolut duct index
    parentDuct->setDuctAbsIndex(nbrDucts);
    nbrDucts++;

    // Set endCAW flag if end duct is reached
    parentDuct->reachedEnd(dLimit);

    // Anchor for transient bronchioles
    if (parentDuct->reachedEnd(dLimit)){

        // Increment number of ending ducts
        nbrEndDucts++;

        // Define shortest duct
        if (lmin > parentDuct->l){
            lmin = parentDuct->l;
            dmin = parentDuct->d;
        }

        return;
    }

    // Grow airway daughters
    parentDuct->grow(0, kappaLb, Phi);
    parentDuct->grow(1, kappaLb, Phi);

    // Call next generation
    genMorphCondAirways(parentDuct->pMajorDaughter);
    genMorphCondAirways(parentDuct->pMinorDaughter);

}


// Compute airway volume
//**************************************************************/
void lung::airwayVolume(duct* parentDuct){

    // Increment total airway volume
    parentDuct->setDuctVolume();
    TAWV += parentDuct->getDuctVolume();

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){

        // Increment total airway length
        Ltot += parentDuct->lc;

        // For lobule template scaling
        totalCompScaling += pow(parentDuct->lc,1.);

        // Set minimum and maximum generation of end ducts
        nbrEndDuctsCheck++;
        if      (nbrEndDuctsCheck == 1)          {endDuctMaxGen = parentDuct->gen;} // descends first to longest
        else if (nbrEndDuctsCheck == nbrEndDucts){endDuctMinGen = parentDuct->gen;}
        return;
    }

    // Call next generation
    airwayVolume(parentDuct->pMajorDaughter);
    airwayVolume(parentDuct->pMinorDaughter);
}


// Size lobules template
//**************************************************************/
void lung::sizeLobuleTemplate(){

    double aAT;

    if (TAWV <= 0.){
        cout << "MISSING VLbLUE: Total airay volume not yet computed.";
        return;
    }

    else if (nbrEndDucts <= 0){
        cout << "MISSING VLbLUE: number of end ducts not yet computed.";
        return;
    }

    else if (totalCompScaling <= 0){
        cout << "MISSING VLbLUE: total compartement scaling by diameter not yet computed.";
        return;
    }

    else if (FRC <= 0.){
        cout << "Functional residual capacity not given.";
        return;
    }

    else if (FRC <= TAWV){
        cout << "NON-PHYSIOLOGIGAL VLbLUE: FRC is smaller than TAWV";
        return;
    }

    else{
        // Size volume for lobules
        if (scalingLbL){
            VLbT = (FRC - TAWV)/(totalCompScaling);
        }
        else{
            VLbT = (FRC - TAWV)/(nbrEndDucts);
        }

        // Lobule volume
        pLbTroot->VLb0 = VLbT;

        // Scaling in duct
        return;
    }
}


// Generate trumpet lobules
//**************************************************************/
void lung::genTrumpetLobuleOnAirways(duct* parentDuct){

    // declarations
    double scaleLbL;

    // Connect lobule to end duct
    if (parentDuct->reachedEnd(dLimit)){

        // Connect a scaled lobule compartement
        parentDuct->connectLobule(scalingLbL, pLbTroot);

        // Set absolute lobule index and increment count for number of ending duct (= number of lobules)
        parentDuct->pLobule->setLobuleAbsIndex(nbrLobules);
        nbrLobules++;

        // Define scale for last duct (used for scaling of trumpe lobule length)
        scaleLbL = (parentDuct->lc*nbrEndDucts)/Ltot;

        // Set trumpet lobule length
        parentDuct->pLobule->setTrumpetLobuleLength(scalingLbL, scaleLbL);

        return;
    }

    // Call next generation
    genTrumpetLobuleOnAirways(parentDuct->pMajorDaughter);
    genTrumpetLobuleOnAirways(parentDuct->pMinorDaughter);

}


// Compute resitance of trumpet lobule
//**************************************************************/
void lung::computeTrumpetLobuleResistance(double d0, double l0, lobule* trumpetLobule){

    // declarations
    double RLb, Ri, li, di;

    // compute lobule resistance as sequence of parallel resistances of equal magnitude
    RLb = 0.;
    for (int k=0; k<maxGenTrumpLob; k++){

        di = d0*pow(kappa,k);
        li = l0*pow(kappa,k);

        Ri = (8.*mu*li)/(pow(di/2,4)*pi);

        RLb += Ri/pow(2.,k);
    }

    // Magnification of total lobular resistance
    // !!! This is a guess and motivated by the fact that the total resistance resulting from a sequence of parallel sub-elements (as computed above) does not inclued any form factors. These may certainly play a role in the in the complex morphology of the human lung.
    RLb *= R_fac;

    trumpetLobule->Racin = RLb;

}


// Generate trumpet lobules
//**************************************************************/
void lung::lobuleVolumeNResistance(duct* parentDuct){

    // declarations
    double d0, l0;

    // correct lobule volume when end duct is reached
    if (parentDuct->reachedEnd(dLimit)){

        // set lobule resistance
        d0 = parentDuct->d;
        l0 = parentDuct->l;
        computeTrumpetLobuleResistance(d0, l0, parentDuct->pLobule);

        // set lobule volume
        //parentDuct->pLobule->setLobuleVolume();
        TLbV += parentDuct->pLobule->getLobuleVolume();


        // calculate stretch distribution in trumpet lobules
        parentDuct->pLobule->calculateStretchDistribution(nbrLobules);

        return;
    }

    // call next generation
    lobuleVolumeNResistance(parentDuct->pMajorDaughter);
    lobuleVolumeNResistance(parentDuct->pMinorDaughter);

}


// Calculate stiffness of Lobule
//**************************************************************/
void lung::calculateStiffnessParametersInTrumpetLobule(duct* parentDuct){


    // correct lobule volume when end duct is reached
    if (parentDuct->reachedEnd(dLimit)){

        // calculate stiffness parameter
        parentDuct->pLobule->calculateStiffnessParameters(nbrLobules, TVm);

        return;
    }

    // call next generation
    calculateStiffnessParametersInTrumpetLobule(parentDuct->pMajorDaughter);
    calculateStiffnessParametersInTrumpetLobule(parentDuct->pMinorDaughter);

}


// Setup gas classes for lung
//**************************************************************/
void lung::setupGaseousSpecies(duct* parentDuct){

    // Declarations
    double hmin, hT, gridsize, alpha, kappa, Gamma, CFL;

    // Grid-size law
    gridsize = (abs(parentDuct->u_max)*dt)/contProp->CFL;

    // Add species (species index, washout index)
    parentDuct->addSpecies(species, washout, gridsize); // (N_2, MBW)

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){

        // Add species in lobule
        parentDuct->pLobule->addSpeciesInTrumpetLobule(species, washout, gridsize);

        // Relate ending duct before returning
        parentDuct->relateSpecies(washout);
        return;
    }

    // Call next generation
    setupGaseousSpecies(parentDuct->pMajorDaughter);
    setupGaseousSpecies(parentDuct->pMinorDaughter);

    // On the way back, relate species parts
    parentDuct->relateSpecies(washout);
}


// Compute time step refinement to satisfy CFL condition in each element
//**************************************************************/
void lung::computeTimeStepRefinement(duct* parentDuct, int opt){

    // declarations
    double Nt_temp, Nt_temp_Lb, u;

    // select advection speed according to 'opt' argument
    if (opt == 0){
        u = abs(parentDuct->u_max);
    }
    if (opt == 1){
        u = abs(parentDuct->u);
    }

    // compute needed timestep refinement in duct
    Nt_temp = (u * dt)/(contProp->CFL * parentDuct->pSpecies0->h);

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){

        // compute needed timestep refinement in lobules (velocity is same as in terminal duct)
        Nt_temp_Lb = (u * dt)/(contProp->CFL * parentDuct->pLobule->pSpecies0->h);

        return;
    }

    // define maximum time step refinement
    if (Nt < Nt_temp){
        Nt = ceil(Nt_temp);

        indAbs_Nt = parentDuct->indAbs;
        gen_Nt = parentDuct->gen;
    }
    if (Nt < Nt_temp_Lb){
        Nt = ceil(Nt_temp_Lb);
    }


    // Call next generation
    computeTimeStepRefinement(parentDuct->pMajorDaughter, opt);
    computeTimeStepRefinement(parentDuct->pMinorDaughter, opt);

}



// Compute total number of grid points
//**************************************************************/
void lung::computeTotalNbrOfGridPoints(duct* parentDuct){

    // Increment count
    NntotD += parentDuct->pSpecies0->N; // nodes
    NntotD++;

    NztotD += parentDuct->pSpecies0->N; // zones

    // Anchor for end duct
    if (parentDuct->reachedEnd(dLimit)){

        // gridpoints in trumpet lobules of longest airway
        NntotLb += parentDuct->pLobule->pSpecies0->N; // nodes
        NntotLb++;

        NztotLb += parentDuct->pLobule->pSpecies0->N; // zones

        return;
    }

    // Call next generation
    computeTotalNbrOfGridPoints(parentDuct->pMajorDaughter);
    computeTotalNbrOfGridPoints(parentDuct->pMinorDaughter);
}


// Initialize VTK arrays for data output
//**************************************************************/
void lung::initializeVTKArrays(){
    /// MESH
    // For ducts
    pZoneTypesD = new int[NztotD];
    pConnectivityD = new int[2*NztotD];
    pCoordinatesD = new float[3*NntotD];

    // For trumpet lobule
    pZoneTypesLb = new int[NztotLb];
    pConnectivityLb = new int[2*NztotLb];
    pCoordinatesLb = new float[3*NntotLb];

    /// OUTPUT VLbRIABLES
    // For ducts
    pNodalLbbsIndD = new float[NntotD];
    pZonalVelocityD = new float[NztotD];
    pNodalPressureD = new float[NntotD];
    pNodalConcentrationID = new float[NntotD];
    pNodalConcentrationIID = new float[NntotD];
    pNodalRadiusD = new float[NntotD];
    isModifiedD = new float[NntotD];

    // For trumpet lobules
    pNodalLbbsIndLb = new float[NntotLb];
    pNodalVelocityLb = new float[NntotLb];
    pNodalConcentrationILb = new float[NntotLb];
    pNodalConcentrationIILb = new float[NntotLb];
    pNodalRadiusLb = new float[NntotLb];
    isModifiedLb = new float[NntotLb];
    pNodalLbcinusVolumeLb = new float[NntotLb];
    pNodalPleuralPressureLb = new float[NntotLb];

    /// OUTPUT VLbRIABLES INFORMATION
    // For duct
    NvarD = 7;

    pDimVarD = new int[NvarD];
    pDimVarD[0] = 1;
    pDimVarD[1] = 1;
    pDimVarD[2] = 1;
    pDimVarD[3] = 1;
    pDimVarD[4] = 1;
    pDimVarD[5] = 1;
    pDimVarD[6] = 1;

    pCenteringD = new int[NvarD];
    pCenteringD[0] = 1;
    pCenteringD[1] = 0;
    pCenteringD[2] = 1;
    pCenteringD[3] = 1;
    pCenteringD[4] = 1;
    pCenteringD[5] = 1;
    pCenteringD[6] = 1;

    pVarnamesD = new const char*[NvarD]; // !!!CAREFULL!!! NOT POSSIBLE PUTTING SPACES
    pVarnamesD[0] = "absind_duct";
    pVarnamesD[1] = "velocity_duct";
    pVarnamesD[2] = "pressure_duct";
    pVarnamesD[3] = "concentrationI_duct";
    pVarnamesD[4] = "concentrationII_duct";
    pVarnamesD[5] = "radius_duct";
    pVarnamesD[6] = "is_modified";

    pAllOutputVariablesD = new float*[NvarD];

    // For trumpet lobule
    NvarLb = 8;

    pDimVarLb = new int[NvarLb];
    pDimVarLb[0] = 1;
    pDimVarLb[1] = 1;
    pDimVarLb[2] = 1;
    pDimVarLb[3] = 1;
    pDimVarLb[4] = 1;
    pDimVarLb[5] = 1;
    pDimVarLb[6] = 1;
    pDimVarLb[7] = 1;

    pCenteringLb = new int[NvarLb];
    pCenteringLb[0] = 1;
    pCenteringLb[1] = 1;
    pCenteringLb[2] = 1;
    pCenteringLb[3] = 1;
    pCenteringLb[4] = 1;
    pCenteringLb[5] = 1;
    pCenteringLb[6] = 1;
    pCenteringLb[7] = 1;

    pVarnamesLb = new const char*[NvarLb]; // !!!CAREFULL!!! NOT POSSIBLE PUTTING SPACES
    pVarnamesLb[0] = "absind_lobules";
    pVarnamesLb[1] = "velocity_lobule";
    pVarnamesLb[2] = "concentrationI_lobule";
    pVarnamesLb[3] = "concentrationII_lobule";
    pVarnamesLb[4] = "radius_lobule";
    pVarnamesLb[5] = "dilatation_lobule";
    pVarnamesLb[6] = "pleural_pressure_lobule";
    pVarnamesLb[7] = "is_modified";

    pAllOutputVariablesLb = new float*[NvarLb];
}


// Initialize matrix for pressure network
//**************************************************************/
void lung::initializeMat(){

    if (nbrDucts <= 0){
        cout << "ERROR: total number of ducts not known yet.";
        return;
    }
    if (bc==0){ // inlet flow boundary conditions => solve for pleural pressure
        coeffMat    = MatrixXd::Zero(nbrDucts+1, nbrDucts+1);
        pressureVec = VectorXd::Zero(nbrDucts+1);
        rhsVec      = VectorXd::Zero(nbrDucts+1);
    }
    else{ // pleural pressure boundary conditions
        coeffMat    = MatrixXd::Zero(nbrDucts, nbrDucts);
        pressureVec = VectorXd::Zero(nbrDucts);
        rhsVec      = VectorXd::Zero(nbrDucts);
    }
    ppl    = -1000;
    dppldt = 0.;

}


// Read the tidal breath and tidal volume data
//**************************************************************/
void lung::readTBTV(breathFlow* allBreathFlow){

    // initialize TB and TV table
    TBTV = MatrixXd::Zero(nbrBreaths,2);

    // read from breath struct
    for (int i = 0; i<nbrBreaths; i++){
        TBTV(i,0) = allBreathFlow[i].TB;
        TBTV(i,1) = allBreathFlow[i].TV;

        // Mean tidal volume and brath period
        TVm += allBreathFlow[i].TV;
        TBm += allBreathFlow[i].TB;
    }

    // Calculate mean
    TVm /= nbrBreaths;
    TBm /= nbrBreaths;
}


// Read Inlet Flow
//**************************************************************/
void lung::readInletFlow(breathFlow* allBreathFlow){

    // Declarartions
    int k;
    double maxInletFlow = 0;

    // Initialize inlet flow vector
    inletFlow = VectorXd::Zero(int(round(Tend/dt)));

    // Read from breath struct
    k = 0;
    for (int i = 0; i<nbrBreaths; i++){
        for (int j = 0; j<allBreathFlow[i].N; j++){
            inletFlow(k) = allBreathFlow[i].flowData[j];

            // Find maximum inlet flow
            if (abs(inletFlow(k)) > maxInletFlow){
                maxInletFlow = abs(inletFlow(k));
                maxIndInletFlow = k;
            }
            k++;
        }
    }

}


// Read Inlet Flow
//**************************************************************/
void lung::readPleuralPressure(breathFlow* allBreathFlow){

    // Declaration
    int k = 0;

    // Initialize inlet flow vector
    pleuralPressure = VectorXd::Zero(int(round(Tend/dt)));

    // Read from breath struct
    k = 0;
    for (int i = 0; i<nbrBreaths; i++){
        for (int j = 0; j<allBreathFlow[i].N; j++){
            pleuralPressure(k) = allBreathFlow[i].pplData[j];
            k++;
        }
    }

    // initialize pleural pressure variable
    ppl = pleuralPressure(0);
}


// Read data for pulsatile correction
//**************************************************************/
void lung::readTransFact(){

    // Check number of lines in in transmissability factor file
    ifstream inFile;
    inFile.open("constant/transFact");
    string unused;
    nbrLinesTransFact = 0;
    while (!inFile.eof()){
        getline(inFile, unused);
        nbrLinesTransFact++;
    }
    inFile.close();
    nbrLinesTransFact--;

    // Initialize transmissibility factor table
    transFact         = MatrixXd::Zero(nbrLinesTransFact, nbrLinesTransFact);
    transFactDomainD  = MatrixXd::Zero(nbrLinesTransFact, nbrLinesTransFact);
    transFactDomainTB = MatrixXd::Zero(nbrLinesTransFact, nbrLinesTransFact);

    // minimum and maximum values of transmissibility domains (diameter, breath period)
    // !!! These values must correspond to the domain range variables hard-coded in the 'pulsatile_transmissibility.py' script used to precompute the transmissibility factor table!!!
    d_min  =  0.0005;
    d_max  =  0.04;
    TB_min =  0.1;
    TB_max = 10.0;

    // Read from text file
    inFile.open("constant/transFact");
    inFile.clear();
    inFile.seekg(0, ios::beg);
    if (inFile.is_open()){
        for (int i = 0; i<nbrLinesTransFact; i++){
            for (int j = 0; j<nbrLinesTransFact; j++){
                inFile >> transFact(i,j);
                transFactDomainTB(i,j) = j*(TB_max - TB_min)/nbrLinesTransFact + TB_min;
                transFactDomainD(i,j)  = i*(d_max  - d_min )/nbrLinesTransFact + d_min;
            }
        }
        inFile.close();
    }
}


// Set tidal volume
//**************************************************************/
void lung::setTBTV(int k){

    // Periode from current breath cylce
    TB = TBTV(k,0);

    // Tidal volume from current breath cylce
    TV = abs(TBTV(k,1));
}


// Read modification table
//**************************************************************/
void lung::readModifications(){

    // check number of lines in in modification file
    ifstream inFile;
    inFile.open("data/modifyLung");
    string unused;
    nbrLinesModTab = 0;
    while (!inFile.eof()){
        getline(inFile, unused);
        nbrLinesModTab++;
    }
    inFile.close();

    // Check for empty last line
    if(unused==""){
        nbrLinesModTab--;
    }

    // Return when modification list is empty
    if (nbrLinesModTab == 0){
        return;
    }

    // Initialize table of modification values
    modTable = MatrixXd::Zero(nbrLinesModTab,5);

    // Read from text file
    inFile.open("data/modifyLung");
    inFile.clear();
    inFile.seekg(0, ios::beg);
    if (inFile.is_open()){
        for (int i = 0; i<nbrLinesModTab; i++){
            inFile >> modTable(i,0); // absolute index (address)
            inFile >> modTable(i,1); // factor duct transmissibility
            inFile >> modTable(i,2); // factor volume modification
            inFile >> modTable(i,3); // factor lobule stretch (affects stiffness)
            inFile >> modTable(i,4); // factor lobule resistance
        }
        inFile.close();
    }

}


// Apply modifications in ducts and lobules
//**************************************************************/
void lung::applyModifications(duct* parentDuct){


    // in conducting airways, apply modifications
    if (nbrLinesModTab > 0){
        // Check weather this duct has to be modified and if yes then do so
        for (int i = 0; i<nbrLinesModTab; i++){
            if (parentDuct->getDuctAbsIndex() == modTable(i,0)){

                // Set duct diameter modification parameter
                parentDuct->fmod *= modTable(i,1);
            }
        }
    }

    // In ending ducts, apply modifications
    if (parentDuct->reachedEnd(dLimit)){

        if (nbrLinesModTab > 0){
            // check whether this duct has to be modified and if yes then do so
            for (int i = 0; i<nbrLinesModTab; i++){
                if (parentDuct->pLobule->getLobuleAbsIndex() == modTable(i,0)){

                    // Modify lobule length
                    parentDuct->pLobule->lLb  *= pow(modTable(i,2), 1./3.);
                    parentDuct->pLobule->lt  *= pow(modTable(i,2), 1./3.);

                    // Modify lobule volume
                    if (int(round(1000*modTable(i,2))) != 1000){
                        // Count number of changed lobules
                        nbrLobulesRed++;
                        TLbVRed += parentDuct->pLobule->VLb*(1-modTable(i,2));
                    }

                    parentDuct->pLobule->VLb  *= modTable(i,2);
                    parentDuct->pLobule->VLb0 *= modTable(i,2);

                    // Modify lobule stretch width
                    parentDuct->pLobule->phi_dV = modTable(i,3);

                    // Modify lobule resistance
                    parentDuct->pLobule->Racin *= modTable(i,4);
                }
            }
        }
        return;
    }

    // Call next generation
    applyModifications(parentDuct->pMajorDaughter);
    applyModifications(parentDuct->pMinorDaughter);

}


// Correct total volume after modification of ducts and lobules
//**************************************************************/
void lung::correctTotalVolumeAfterModification(duct* parentDuct){

    // Declarations
    bool wasModified;

    // In ending duct apply correction
    if (parentDuct->reachedEnd(dLimit)){

        // Check whether the lobule volume has been modifed (then must not be corrected)
        if (nbrLinesModTab > 0){

            for (int i = 0; i<nbrLinesModTab; i++){
                if (parentDuct->pLobule->getLobuleAbsIndex() == modTable(i,0)){


                    if (int(round(1000*modTable(i,2))) != 1000){
                        wasModified = true;
                    }
                    else {
                        wasModified = false;
                    }
                }
            }
        }
        else{
            wasModified = false;
        }

        // Correct the remaining lobules
        if (!wasModified){
            parentDuct->pLobule->VLb  += TLbVRed/(nbrLobules - nbrLobulesRed);
            parentDuct->pLobule->VLb0 += TLbVRed/(nbrLobules - nbrLobulesRed);
        }
        TLbV_aft_mod += parentDuct->pLobule->VLb;
        return;
    }

    // Call next generation
    correctTotalVolumeAfterModification(parentDuct->pMajorDaughter);
    correctTotalVolumeAfterModification(parentDuct->pMinorDaughter);
}


// Set transmissibility in ducts
//**************************************************************/
void lung::airwaySetTK1(duct* parentDuct){

    // declarations
    double TF;
    double evalTB, evalD;

    // Set transmissibility in current duct
    evalTB = TB;
    evalD  = parentDuct->d;
    TF = parentDuct->bilinearInterpolation(evalTB, evalD, nbrLinesTransFact, transFact, transFactDomainTB, transFactDomainD);
    parentDuct->setTransmissibility(mu, TF);

    // Also set K1 in end ducts
    if (parentDuct->reachedEnd(dLimit)){
        parentDuct->setK1(dt/Nt);
        return;
    }

    // Call next generation
    airwaySetTK1(parentDuct->pMajorDaughter);
    airwaySetTK1(parentDuct->pMinorDaughter);
}


// Write coefficient matrix
//**************************************************************/
void lung::writeCoeffMatrix(duct* parentDuct){

    // Write entries for equation number indAbs
    parentDuct->writeCoeff(bc, nbrDucts, coeffMat);

    if (parentDuct->reachedEnd(dLimit)){
        return;
    }

    // Call next generation
    writeCoeffMatrix(parentDuct->pMajorDaughter);
    writeCoeffMatrix(parentDuct->pMinorDaughter);
}


// Interpolate input data (inlet flow rate or pleureal pressure) for sub-time step
//**************************************************************/
void lung::interpInputData(int it, int n){

    // Declarations
    double Qin0, Qinp1, ppl0, pplp1;

    if (bc==0){
        // Interpolate inlet flow value for sub-time step
        Qin0 = inletFlow(it);
        if (it == int(round(Tend/dt)) - 1){
            Qinp1 = 0;
        }
        else{
            Qinp1 = inletFlow(it+1);
        }

        Qin = Qin0 + double(n)/Nt*(Qinp1 - Qin0);
    }
    else{
        // Interpolate pleural pressure value for sub-time step
        ppl0 = pleuralPressure(it);
        if (it == 0){
            pplp1  = pleuralPressure(it+1);
            dppldt = (pplp1 - ppl0)/dt;
        }
        else if (it == int(round(Tend/dt)) - 1){
            pplp1  = pleuralPressure(0);
            dppldt = (pplp1 - ppl0)/dt;
        }
        else{
            pplp1  = pleuralPressure(it+1);
            dppldt = (pplp1 - pleuralPressure(it-1))/(2.*dt);
        }

        ppl = ppl0 + double(n)/Nt*(pplp1 - ppl0);
    }
}


// Write right-hand-side for pressure network
//**************************************************************/
void lung::writeRHSVector(duct* parentDuct, int opt){

    // declarations
    double Qin_, varinp;

    if (opt==0){ // regular case during simulation run
        if (bc == 0){
            Qin_ = Qin;
            varinp = ppl;
        }
        else{
            Qin_ = 0.; // not needed when using pleural pressure boundary conditions
            varinp = dppldt;
        }
    }
    else{ // prior to simulation start for flow distr. estimation
        if (bc==0){
            Qin_ = inletFlow(maxIndInletFlow);
            varinp = -2500.; // as minimum pleural pressure
        }
        else{
            Qin_ = TVm*pi/TBm;
            varinp = -2500.; // as minimum pleural pressure
        }
    }

    // Write entries for equation in number indAbs
    parentDuct->writeRHS(bc, nbrEndDucts, TVm, dt/Nt, patm, varinp, Qin_, rhsVec);

    if (parentDuct->reachedEnd(dLimit)){
        return;
    }

    // Call next generation
    writeRHSVector(parentDuct->pMajorDaughter, opt);
    writeRHSVector(parentDuct->pMinorDaughter, opt);
}


// Distribute pressure over network
//**************************************************************/
void lung::distrPressure(duct* parentDuct){

    // Set pressure
    parentDuct->p = pressureVec(parentDuct->indAbs);

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){
        return;
    }

    // Call next generation
    distrPressure(parentDuct->pMajorDaughter);
    distrPressure(parentDuct->pMinorDaughter);
}


// Compute Flow
//**************************************************************/
void lung::computeFlow(duct* parentDuct){

    // Set flow
    parentDuct->setFlow(washout, patm);

    if (parentDuct->gen==0 && bc==1){
        Qin = parentDuct->Q;
    }

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){

        // transfer to inlet flow / velocity in lobule
        parentDuct->pLobule->Q = parentDuct->Q;
        parentDuct->pLobule->u = parentDuct->u;

        // increment total outlet flow
        Qouttot += parentDuct->Q;
        return;
    }

    // call next generation
    computeFlow(parentDuct->pMajorDaughter);
    computeFlow(parentDuct->pMinorDaughter);
}


// Estimate flow for CFL Condition
//**************************************************************/
void lung::estimateFlow(duct* parentDuct){

    // Estimate flow
    parentDuct->estimateFlow(washout, patm, Qin);

    // Set smallest absolute velocity in morphology
    if (parentDuct->l == lmin){
        umin = abs(parentDuct->u);
    }

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){
        return;
    }

    // Call next generation
    estimateFlow(parentDuct->pMajorDaughter);
    estimateFlow(parentDuct->pMinorDaughter);
}


// Update flow properties in lobule
//**************************************************************/
void lung::updateLobule(duct* parentDuct){

    // Declarations
    double Flow, volumeRatio;

    // Update lobule volume when end duct is reached
    if (parentDuct->reachedEnd(dLimit)){

        // Set and get ratio of updated volume to former volume
        Flow = parentDuct->Q;
        parentDuct->pLobule->setVolumeRatio(dt/Nt, Flow);
        volumeRatio = parentDuct->pLobule->getVolumeRatio();

        // Update inlet flow and volume of trumpet lobule
        parentDuct->pLobule->updateTrumpetLobule(washout, dt/Nt, Flow);
        return;
    }

    // Call next generation
    updateLobule(parentDuct->pMajorDaughter);
    updateLobule(parentDuct->pMinorDaughter);
}


// Time integration of concentration in ducts and lobules
//**************************************************************/
void lung::updateConcentrationInDucts(duct* parentDuct){

    // Declarations
    bool iT, rE;
    double cCondAW0, cCondAW1, c1Lobule0, c2Lobule0, c1Lobule1, c2Lobule1, dcdxCondAW0, dcdxCondAW1, hCondAW0, hCondAW1;

    // Pick arguments for different cases
    if (parentDuct->gen==0){
        iT = true;
        rE = false;
        c1Lobule0 = c2Lobule0 = c1Lobule1 = c2Lobule1 = -1.; // not used
    }
    else if (parentDuct->reachedEnd(dLimit)){
        iT = false;
        rE = true;
        if (washout==0){
            c1Lobule0 = parentDuct->pLobule->pSpecies0->c1A;
            c2Lobule0 = parentDuct->pLobule->pSpecies0->c2A;
        }
        if (washout==1){
            c1Lobule0 = parentDuct->pLobule->pSpecies0->c1A;
            c2Lobule0 = parentDuct->pLobule->pSpecies0->c2A;
            c1Lobule1 = parentDuct->pLobule->pSpecies1->c1A;
            c2Lobule1 = parentDuct->pLobule->pSpecies1->c2A;
        }
    }
    else{
        iT = false;
        rE = false;
        c1Lobule0 = c2Lobule0 = c1Lobule1 = c2Lobule1 = -1.; // not used
    }

    // Update concentration in current duct
    if (washout==0){
        parentDuct->pSpecies0->updateConcentrationInDuct(iT, rE, dt/Nt, cin, c1Lobule0, c2Lobule0);
    }
    if (washout==1){
        parentDuct->pSpecies0->updateConcentrationInDuct(iT, rE, dt/Nt, cin, c1Lobule0, c2Lobule0);
        parentDuct->pSpecies1->updateConcentrationInDuct(iT, rE, dt/Nt, cin, c1Lobule1, c2Lobule1);
    }

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){

        // Define values to be passed to lobule
        if (washout==0){
            cCondAW0 = parentDuct->pSpecies0->cCAW;
            dcdxCondAW0 = parentDuct->pSpecies0->dcdxCAW;
            hCondAW0 = parentDuct->pSpecies0->h;
        }
        if (washout==1){
            cCondAW0 = parentDuct->pSpecies0->cCAW;
            dcdxCondAW0 = parentDuct->pSpecies0->dcdxCAW;
            hCondAW0 = parentDuct->pSpecies0->h;

            cCondAW1 = parentDuct->pSpecies1->cCAW;
            dcdxCondAW1 = parentDuct->pSpecies1->dcdxCAW;
            hCondAW1 = parentDuct->pSpecies1->h;
        }

        // Update concentration in trumpet lobule
        if (washout==0){
            parentDuct->pLobule->pSpecies0->updateConcentrationInTrumpetLobule(dt/Nt, cCondAW0, dcdxCondAW0, hCondAW0);
        }
        if (washout==1){
            parentDuct->pLobule->pSpecies0->updateConcentrationInTrumpetLobule(dt/Nt, cCondAW0, dcdxCondAW0, hCondAW0);
            parentDuct->pLobule->pSpecies1->updateConcentrationInTrumpetLobule(dt/Nt, cCondAW1, dcdxCondAW1, hCondAW1);
        }
        return;
    }

    // Call next generation
    updateConcentrationInDucts(parentDuct->pMajorDaughter);
    updateConcentrationInDucts(parentDuct->pMinorDaughter);

}


// Collect data in trumpet lobule
//**************************************************************/
void lung::collectAllOutputDataInTrumpetLobule(duct* parentDuct){

    // Declarations
    int N, fni, fzi, fcri, fcni, k2, k3;
    double x, y, z, dx, dy, dz, dp;
    double dist, lLb, L, a;

    N = parentDuct->pLobule->pSpecies0->N;
    fni = firstNodeIndLb;
    fzi = firstZoneIndLb;
    fcri = firstCoordIndLb;
    fcni = firstConnIndLb;

    // Node coordinates
    for (int k=0; k<=N; k++){
        k3 = 3*k;

        lLb = parentDuct->pLobule->lLb;
        L  = parentDuct->pLobule->L;

        dist = k*lLb/double(N) / 5.;
        x = parentDuct->x[1] + cos(parentDuct->phi_z)*cos(parentDuct->phi)*dist;
        y = parentDuct->y[1] + cos(parentDuct->phi_z)*sin(parentDuct->phi)*dist;
        z = parentDuct->z[1] + sin(parentDuct->phi_z)*dist;


        pCoordinatesLb[fcri+k3]   = x;
        pCoordinatesLb[fcri+k3+1] = y;
        pCoordinatesLb[fcri+k3+2] = z;
    }

    // Zone types
    for (int k=0; k<N; k++){
        pZoneTypesLb[fzi+k] = 3; // is VISIT_LINE (see visit_writer.h)
    }

    // Connectivity
    for (int k=0; k<N; k++){
        k2 = 2*k;

        pConnectivityLb[fcni+k2] = fni + k;
        pConnectivityLb[fcni+k2+1] = fni + k + 1;
    }


    // Nodal data
    for (int k=0; k<=N; k++){

        // Absolut index
        pNodalLbbsIndLb[fni+k] = parentDuct->pLobule->indAbsLob;

        // Concentration values
        if (washout==0){
            pNodalConcentrationILb[fni+k]  = parentDuct->pLobule->pSpecies0->C(k);
            pNodalConcentrationIILb[fni+k] = -1;
        }
        if (washout==1){
            pNodalConcentrationILb[fni+k]  = parentDuct->pLobule->pSpecies0->C(k);
            pNodalConcentrationIILb[fni+k] = parentDuct->pLobule->pSpecies1->C(k);
        }

        if (parentDuct->pLobule->pSpecies0->C(k) > 1.01){
          cout << "ERROR: Concentration in trumpet lobule higher than 1: C(k) = " << parentDuct->pLobule->pSpecies0->C(k) << endl;
        }

        // Radius (scaling factor 0.1 is for better visibility when displayed in a visualization tool)
        pNodalRadiusLb[fni+k] = 0.1*sqrt(4. * parentDuct->pLobule->pSpecies0->ATr(k)/pi)/2.;

        // Velocity
        pNodalVelocityLb[fni+k] = parentDuct->pLobule->pSpecies0->UphTr(k);

        // Lobule volume
        pNodalLbcinusVolumeLb[fni+k] = parentDuct->pLobule->VLb - parentDuct->pLobule->VLb0;

        // Pleural pressure
        pNodalPleuralPressureLb[fni+k] = ppl;
    }


    // Increment indices;
    firstNodeIndLb += (N + 1);
    firstZoneIndLb +=  N;
    firstCoordIndLb += 3*(N + 1);
    firstConnIndLb += 2*N;

    return;

}


// Collect all output data
//**************************************************************/
void lung::collectAllOutputData(duct* parentDuct){

    // Declarations
    int N, fni, fzi, fcri, fcni, k2, k3;
    double x, y, z, dx, dy, dz, dp;

    N = parentDuct->pSpecies0->N;
    fni = firstNodeIndD;
    fzi = firstZoneIndD;
    fcri = firstCoordIndD;
    fcni = firstConnIndD;

    // Node coordinates
    for (int k=0; k<=N; k++){
        k3 = 3*k;

        dx = parentDuct->x[1] - parentDuct->x[0];
        dy = parentDuct->y[1] - parentDuct->y[0];
        dz = parentDuct->z[1] - parentDuct->z[0];

        x = parentDuct->x[0] + k*dx/double(N);
        y = parentDuct->y[0] + k*dy/double(N);
        z = parentDuct->z[0] + k*dz/double(N);

        pCoordinatesD[fcri+k3]   = x;
        pCoordinatesD[fcri+k3+1] = y;
        pCoordinatesD[fcri+k3+2] = z;
    }

    // Zone types
    for (int k=0; k<N; k++){
        pZoneTypesD[fzi+k] = 3; // is VISIT_LINE (see visit_writer.h)
    }

    // Connectivity
    for (int k=0; k<N; k++){
        k2 = 2*k;

        pConnectivityD[fcni+k2] = fni + k;
        pConnectivityD[fcni+k2+1] = fni + k + 1;
    }


    // Nodal data
    for (int k=0; k<=N; k++){

        // Absolut index
        pNodalLbbsIndD[fni+k] = parentDuct->indAbs;

        // Pressure
        if (parentDuct->gen == 0) { // trachea
            dp = parentDuct->p - patm;
            pNodalPressureD[fni+k] = patm + k*dp/double(N);
        }
        else{ // inside lung morphology
            dp = parentDuct->p - parentDuct->pParent->p;
            pNodalPressureD[fni+k] = parentDuct->pParent->p + k*dp/double(N);
        }

        // Concentration values
        if (washout==0){
            pNodalConcentrationID[fni+k]  = parentDuct->pSpecies0->C(k);
            pNodalConcentrationIID[fni+k] = -1;
        }
        if (washout==1){
            pNodalConcentrationID[fni+k]  = parentDuct->pSpecies0->C(k);
            pNodalConcentrationIID[fni+k] = parentDuct->pSpecies1->C(k);
        }

        if (parentDuct->pSpecies0->C(k) > 1.01){
          cout << "ERROR: Concentration in duct higher than 1: C(k) = " << parentDuct->pSpecies0->C(k) << endl;
        }

        // Radius (scaling factor 0.1 is for better visibility when displayed in a visualization tool)
        pNodalRadiusD[fni+k] = 0.1*parentDuct->d/2.;


    }


    // Zonal data
    for (int k=0; k<N; k++){

        // Velocity
        pZonalVelocityD[fzi+k] = parentDuct->u;
    }


    // Increment indices;
    firstNodeIndD += (N + 1);
    firstZoneIndD +=  N;
    firstCoordIndD += 3*(N + 1);
    firstConnIndD += 2*N;

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){

        // Collect in trumpet lobule
        collectAllOutputDataInTrumpetLobule(parentDuct);

        return;
    }

    // Call next generation
    collectAllOutputData(parentDuct->pMajorDaughter);
    collectAllOutputData(parentDuct->pMinorDaughter);

}


// Write full lung data
//**************************************************************/
void lung::writeFullLungData(duct* parentDuct){

    // TOTAL LUNG OUTPUTS (VTK)
    // Collect data
    firstNodeIndD = firstZoneIndD = firstCoordIndD = firstConnIndD = 0;
    firstNodeIndLb = firstZoneIndLb = firstCoordIndLb = firstConnIndLb = 0;
    collectAllOutputData(parentDuct);

    // Write duct data
    pAllOutputVariablesD[0] = pNodalLbbsIndD;
    pAllOutputVariablesD[1] = pZonalVelocityD;
    pAllOutputVariablesD[2] = pNodalPressureD;
    pAllOutputVariablesD[3] = pNodalConcentrationID;
    pAllOutputVariablesD[4] = pNodalConcentrationIID;
    pAllOutputVariablesD[5] = pNodalRadiusD;
    pAllOutputVariablesD[6] = isModifiedD;


    char filenameD[100];
    sprintf(filenameD, "data/duct/TotalLungOutputDuct%05d.vtk", frame);

    write_unstructured_mesh(filenameD, 0, NntotD, pCoordinatesD, NztotD, pZoneTypesD, pConnectivityD, NvarD, pDimVarD, pCenteringD, pVarnamesD, pAllOutputVariablesD);

    // Write lobule data
    pAllOutputVariablesLb[0] = pNodalLbbsIndLb;
    pAllOutputVariablesLb[1] = pNodalVelocityLb;
    pAllOutputVariablesLb[2] = pNodalConcentrationILb;
    pAllOutputVariablesLb[3] = pNodalConcentrationIILb;
    pAllOutputVariablesLb[4] = pNodalRadiusLb;
    pAllOutputVariablesLb[5] = pNodalLbcinusVolumeLb;
    pAllOutputVariablesLb[6] = pNodalPleuralPressureLb;
    pAllOutputVariablesLb[7] = isModifiedLb;

    char filenameLb[100];
    sprintf(filenameLb, "data/lobule/TotalLungOutputLobule%05d.vtk", frame);

    write_unstructured_mesh(filenameLb, 0, NntotLb, pCoordinatesLb, NztotLb, pZoneTypesLb, pConnectivityLb, NvarLb, pDimVarLb, pCenteringLb, pVarnamesLb, pAllOutputVariablesLb);
    frame += 1;
}


// Testing issues
//**************************************************************/
void lung::writeOut(duct* parentDuct){

    // Declarations
    int gen;
    int fb_case;

    double Q_p, Q_maj, Q_min;
    double u_p, u_maj, u_min;
    double d_p, d_maj, d_min;
    double flow_balance, flow_balance_rel;
    double fb_p, fb_maj, fb_min;

    MatrixXd lo, d1, d2;
    VectorXd i, rhs, b;

    // Current generation
    gen = parentDuct->gen;

    if (!parentDuct->reachedEnd(dLimit)){
        u_p = parentDuct->pSpecies0->u;
        u_maj = parentDuct->pMajorDaughter->pSpecies0->u;
        u_min = parentDuct->pMinorDaughter->pSpecies0->u;

        d_p   = sqrt(4*parentDuct->pSpecies0->S/pi);
        d_min = sqrt(4*parentDuct->pMinorDaughter->pSpecies0->S/pi);
        d_maj = sqrt(4*parentDuct->pMajorDaughter->pSpecies0->S/pi);

        fb_p = parentDuct->pSpecies0->fb_p;
        fb_maj = parentDuct->pMajorDaughter->pSpecies0->fb_maj;
        fb_min = parentDuct->pMinorDaughter->pSpecies0->fb_min;

        // Flux bilance
        Q_p = parentDuct->pSpecies0->Q;
        Q_maj = parentDuct->pMajorDaughter->pSpecies0->Q;
        Q_min = parentDuct->pMinorDaughter->pSpecies0->Q;

        flow_balance = Q_p - Q_maj - Q_min;
        flow_balance_rel = flow_balance/Q_p;

        if (Q_p >= 0 && Q_maj >= 0 && Q_min >= 0){
            fb_case = 0;
        }
        else if (Q_p >= 0 && Q_maj < 0 && Q_min >= 0){
            fb_case = 1;
        }
        else if (Q_p >= 0 && Q_maj >= 0 && Q_min < 0){
            fb_case = 2;
        }
        else if (Q_p < 0 && Q_maj >= 0 && Q_min < 0){
            fb_case = 3;
        }
        else if (Q_p < 0 && Q_maj < 0 && Q_min >= 0){
            fb_case = 4;
        }
        else if (Q_p < 0 && Q_maj < 0 && Q_min < 0){
            fb_case = 5;
        }
        else {
            fb_case = 6;
        }

    }

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){
        return;
    }

    // Call next generation
    writeOut(parentDuct->pMajorDaughter);
}


// Write primary lung data
void lung::writePrimaryLungData(int it, breathFlow* allBreathFlow, breathResults* allBreathResults, duct* parentDuct){

    // Declaration
    int breath_ind, data_ind;
    double TB_past;

    // Calculate Index
    breath_ind = TB_past = 0;
    for (int m=0; m<=it; m++){

        if (m == int(round(TB_past/dt))+int(round(TBTV(breath_ind,0)/dt))) {
            TB_past += TBTV(breath_ind,0);
            breath_ind++;
        }
    }

    // index for data point in current time step
    data_ind = it - int(round(TB_past/dt));

    // write time
    allBreathResults[breath_ind].time[data_ind] = it*dt;
    if (it > 6470 && it <= 6480){
        //cout << "it = " << it << "; dt = " << dt << "; data_ind = " << data_ind  << "; TB_past = " << TB_past << "; int(round(TB_past/dt)) = " << int(round(TB_past/dt)) << endl;
    }

    // check whether concentration is bounded
    if (parentDuct->pSpecies0->C(0) > 1.01){
      cout << "ERROR: Concentration at inlet higher than 1: C_in = " << parentDuct->pSpecies0->C(0) << endl;
    }

    // Write species concentration at inlet
    if (washout==0){
        allBreathResults[breath_ind].speciesIData[data_ind] = parentDuct->pSpecies0->C(0);
        allBreathResults[breath_ind].speciesIIData[data_ind] = -1;
    }
    if (washout==1){
        allBreathResults[breath_ind].speciesIData[data_ind]  = parentDuct->pSpecies0->C(0);
        allBreathResults[breath_ind].speciesIIData[data_ind] = parentDuct->pSpecies1->C(0);
    }

    // write pleural pressure
    allBreathResults[breath_ind].pleuralPressure[data_ind] = ppl;

}


// Update old concentration vector in duct
//**************************************************************/
void lung::updateOldConc(duct* parentDuct){

    parentDuct->pSpecies0->updateOldConc();

        if (parentDuct->reachedEnd(dLimit)){
            return;
        }

    updateOldConc(parentDuct->pMajorDaughter);
    updateOldConc(parentDuct->pMinorDaughter);

}

#include "lung.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <eigen_3_3_4/Dense>



// Constructor
//**************************************************************/
lung::lung(init *initData, controlProperties *contProp_,
            systemProperties *sysProp_, transportProperties *transProp_){

    pi = atan(1.0)*4.0;

    // Write class members
    //**********************************************************/
    contProp            = contProp_;
    sysProp             = sysProp_;
    transProp           = transProp_;


    // Initialize lung variables and pass values
    //**********************************************************/
    dLimit              = initData->dL;
    FRC                 = initData->FRC;
    washout             = initData->washout;
    species             = initData->species;
    dt                  = initData->dt;
    Tend                = initData->Tend;
    nbrBreaths          = initData->nbrBreaths;


    // Set bifurcation parameters
    maxGenTrumpAcin     = sysProp->maxGenAc;
    patm                = sysProp->p0;
    r                   = sysProp->r;
    eta                 = sysProp->eta;
    kappaA[0]           = pow((1.-r),(1./eta));
    kappaA[1]           = pow((r),(1./eta));

    // Set homothety ratio in the trumpet lobule (kappa_hat for cross-section)
    kappa               = sysProp->kappa;
    kappa_hat          = 2.*kappa*kappa;

    // Set inlet concentration
    cin                 = sysProp->cin;
    bc                  = sysProp->bc;

    // Set transport properties
    mu                  = transProp->mu;
    rho                 = transProp->rho;

    // magnification factor for total lobular resistance
    R_fac               = sysProp->R_fac;

    // flag for trumpet lobule scaling
    scalingTL           = bool(sysProp->scalingTL);


    // Initialize remaining variables
    //**********************************************************/
    Nt = 1;
    indAbs_Nt = gen_Nt = 0;

    TVm = 0;
    TBm = 0;
    E          = 1e09;

    NntotD = NztotD = NntotA = NztotA = 0;

    firstNodeIndD = firstZoneIndD = firstCoordIndD = firstConnIndD = 0;
    firstNodeIndA = firstZoneIndA = firstCoordIndA = firstConnIndA = 0;

    umin       = 0;

    // Set input flow
    Qin        = 0.;
    Qouttot    = 0.;

    TAWV = TAV = TAV_aft_mod = TRLV = TAVRed = 0;
    totalCompScaling = 0.;
    nbrEndDucts = nbrDucts = nbrEndDuctsCheck = endDuctMinGen = endDuctMinGen = 0;
    nbrAcini = nbrAciniRed = 0;

    // Initialize Trachea
    //**********************************************************/
    pTrachea = new duct(contProp, sysProp, transProp);

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

    // Initialize acinus template
    //**********************************************************/
    pATroot = new acinus(contProp, sysProp, transProp);

    // Acinus volume of template
    VAAT = 0;

    // Compliance of acini
    pATroot->E = E;

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
    delete[] pZoneTypesA;
    delete[] pConnectivityD;
    delete[] pConnectivityA;

    delete[] pCoordinatesD;
    delete[] pCoordinatesA;

    delete[] pZonalVelocityD;
    delete[] pNodalVelocityA;
    delete[] pNodalPressureD;
    delete[] pNodalConcentrationID;
    delete[] pNodalConcentrationIA;
    delete[] pNodalConcentrationIID;
    delete[] pNodalConcentrationIIA;

    delete[] pDimVarD;
    delete[] pDimVarA;
    delete[] pCenteringD;
    delete[] pCenteringA;
    delete[] pVarnamesD;
    delete[] pVarnamesA;

    delete[] pAllOutputVariablesD;
    delete[] pAllOutputVariablesA;

    // Create rubber that moves through the duct tree and deletes its branches
    duct* rubberD = new duct(contProp, sysProp, transProp);

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

        // Create rubber that moves through the acinus tree and delets its branches
        acinus* rubberA = new acinus(contProp, sysProp, transProp);

        // Delete trumpet acinus compartement
        rubberA = parentDuct->pAcinus;

        delete rubberA;

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
    parentDuct->grow(0, kappaA, Phi);
    parentDuct->grow(1, kappaA, Phi);

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

        // For acinus template scaling
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


// Size acini template
//**************************************************************/
void lung::sizeAcinusTemplate(){

    double aAT;

    if (TAWV <= 0.){
        cout << "MISSING VALUE: Total airay volume not yet computed.";
        return;
    }

    else if (nbrEndDucts <= 0){
        cout << "MISSING VALUE: number of end ducts not yet computed.";
        return;
    }

    else if (totalCompScaling <= 0){
        cout << "MISSING VALUE: total compartement scaling by diameter not yet computed.";
        return;
    }

    else if (FRC <= 0.){
        cout << "Functional residual capacity not given.";
        return;
    }

    else if (FRC <= TAWV){
        cout << "NON-PHYSIOLOGIGAL VALUE: FRC is smaller than TAWV";
        return;
    }

    else{
        // Size volume for acini
        if (scalingTL){
            VAAT = (FRC - TAWV)/(totalCompScaling);
        }
        else{
            VAAT = (FRC - TAWV)/(nbrEndDucts);
        }

        // Acinus volume
        pATroot->VA0 = VAAT;

        // Scaling in duct
        return;
    }
}


// Generate trumpet acini
//**************************************************************/
void lung::genTrumpetAcinusOnAirways(duct* parentDuct){

    // declarations
    double scaleTL;

    // Connect acinus to end duct
    if (parentDuct->reachedEnd(dLimit)){

        // Connect a scaled acinus compartement
        parentDuct->connectAcinus(scalingTL, pATroot);

        // Set absolute acinus index and increment count for number of ending duct (= number of acini)
        parentDuct->pAcinus->setAcinusAbsIndex(nbrAcini);
        nbrAcini++;

        // Define scale for last duct (used for scaling of trumpe lobule length)
        scaleTL = (parentDuct->lc*nbrEndDucts)/Ltot;

        // Set trumpet acinus length
        parentDuct->pAcinus->setTrumpetLength(scalingTL, scaleTL);

        return;
    }

    // Call next generation
    genTrumpetAcinusOnAirways(parentDuct->pMajorDaughter);
    genTrumpetAcinusOnAirways(parentDuct->pMinorDaughter);

}


// Compute resitance of trumpet acinus
//**************************************************************/
void lung::computeTrumpetAcinusResistance(double d0, double l0, acinus* trumpetAcinus){

    // declarations
    double RA, Ri, li, di;

    // compute acinus resistance as sequence of parallel resistances of equal magnitude
    RA = 0.;
    for (int k=0; k<maxGenTrumpAcin; k++){

        di = d0*pow(kappa,k);
        li = l0*pow(kappa,k);

        Ri = (8.*mu*li)/(pow(di/2,4)*pi);

        RA += Ri/pow(2.,k);
    }

    // Magnification of total lobular resistance
    // !!! This is a pure guess and motivated by the fact that the total resistance resulting from a sequence of parallel sub-elements (as computed above) does not inclued any form factors. These may certainly play a role in the in the complex morphology of the human lung.
    RA *= R_fac;

    trumpetAcinus->Racin = RA;

}


// Generate trumpet acini
//**************************************************************/
void lung::acinusVolumeNResistance(duct* parentDuct){

    // declarations
    double d0, l0;

    // correct acinus volume when end duct is reached
    if (parentDuct->reachedEnd(dLimit)){

        // set acinus resistance
        d0 = parentDuct->d;
        l0 = parentDuct->l;
        computeTrumpetAcinusResistance(d0, l0, parentDuct->pAcinus);

        // set acinus volume
        //parentDuct->pAcinus->setAcinusVolume();
        TAV += parentDuct->pAcinus->getAcinusVolume();


        // calculate stretch distribution in trumpet acinus
        parentDuct->pAcinus->calculateStretchDistribution(nbrAcini);

        return;
    }

    // call next generation
    acinusVolumeNResistance(parentDuct->pMajorDaughter);
    acinusVolumeNResistance(parentDuct->pMinorDaughter);

}


// Calculate stiffness of Acinus
//**************************************************************/
void lung::calculateStiffnessParametersInTrumpetAcinus(duct* parentDuct){


    // correct acinus volume when end duct is reached
    if (parentDuct->reachedEnd(dLimit)){

        // calculate stiffness parameter
        parentDuct->pAcinus->calculateStiffnessParameters(nbrAcini, TVm);

        return;
    }

    // call next generation
    calculateStiffnessParametersInTrumpetAcinus(parentDuct->pMajorDaughter);
    calculateStiffnessParametersInTrumpetAcinus(parentDuct->pMinorDaughter);

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

        // Add species in acinus
        parentDuct->pAcinus->addSpeciesInTrumpet(species, washout, gridsize);

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
    double Nt_temp, Nt_temp_A, u;

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

        // compute needed timestep refinement in acinus (velocity is same as in terminal duct)
        Nt_temp_A = (u * dt)/(contProp->CFL * parentDuct->pAcinus->pSpecies0->h);

        return;
    }

    // define maximum time step refinement
    if (Nt < Nt_temp){
        Nt = ceil(Nt_temp);

        indAbs_Nt = parentDuct->indAbs;
        gen_Nt = parentDuct->gen;
    }
    if (Nt < Nt_temp_A){
        Nt = ceil(Nt_temp_A);
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

        // gridpoints in trumpet acinus of longest airway
        NntotA += parentDuct->pAcinus->pSpecies0->N; // nodes
        NntotA++;

        NztotA += parentDuct->pAcinus->pSpecies0->N; // zones

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

    // For trumpet acinus
    pZoneTypesA = new int[NztotA];
    pConnectivityA = new int[2*NztotA];
    pCoordinatesA = new float[3*NntotA];

    /// OUTPUT VARIABLES
    // For ducts
    pNodalAbsIndD = new float[NntotD];
    pZonalVelocityD = new float[NztotD];
    pNodalPressureD = new float[NntotD];
    pNodalConcentrationID = new float[NntotD];
    pNodalConcentrationIID = new float[NntotD];
    pNodalRadiusD = new float[NntotD];
    isModifiedD = new float[NntotD];

    // For trumpet acinus
    pNodalAbsIndA = new float[NntotA];
    pNodalVelocityA = new float[NntotA];
    pNodalConcentrationIA = new float[NntotA];
    pNodalConcentrationIIA = new float[NntotA];
    pNodalRadiusA = new float[NntotA];
    isModifiedA = new float[NntotA];
    pNodalAcinusVolumeA = new float[NntotA];
    pNodalPleuralPressureA = new float[NntotA];

    /// OUTPUT VARIABLES INFORMATION
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

    // For trumpet acinus
    NvarA = 8;

    pDimVarA = new int[NvarA];
    pDimVarA[0] = 1;
    pDimVarA[1] = 1;
    pDimVarA[2] = 1;
    pDimVarA[3] = 1;
    pDimVarA[4] = 1;
    pDimVarA[5] = 1;
    pDimVarA[6] = 1;
    pDimVarA[7] = 1;

    pCenteringA = new int[NvarA];
    pCenteringA[0] = 1;
    pCenteringA[1] = 1;
    pCenteringA[2] = 1;
    pCenteringA[3] = 1;
    pCenteringA[4] = 1;
    pCenteringA[5] = 1;
    pCenteringA[6] = 1;
    pCenteringA[7] = 1;

    pVarnamesA = new const char*[NvarA]; // !!!CAREFULL!!! NOT POSSIBLE PUTTING SPACES
    pVarnamesA[0] = "absind_acinus";
    pVarnamesA[1] = "velocity_acinus";
    pVarnamesA[2] = "concentrationI_acinus";
    pVarnamesA[3] = "concentrationII_acinus";
    pVarnamesA[4] = "radius_acinus";
    pVarnamesA[5] = "dilatation_acinus";
    pVarnamesA[6] = "pleural_pressure_acinus";
    pVarnamesA[7] = "is_modified";

    pAllOutputVariablesA = new float*[NvarA];
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
    inFile.open("data/transFact");
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
    inFile.open("data/transFact");
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
            inFile >> modTable(i,3); // factor acinus stretch (affects stiffness)
            inFile >> modTable(i,4); // factor acinus resistance
        }
        inFile.close();
    }

}


// Apply modifications in ducts and acini
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
                if (parentDuct->pAcinus->getAcinusAbsIndex() == modTable(i,0)){

                    // Modify acinus length
                    parentDuct->pAcinus->lA  *= pow(modTable(i,2), 1./3.);
                    parentDuct->pAcinus->lt  *= pow(modTable(i,2), 1./3.);

                    // Modify acinus volume
                    if (int(round(1000*modTable(i,2))) != 1000){
                        // Count number of changed acini
                        nbrAciniRed++;
                        TAVRed += parentDuct->pAcinus->VA*(1-modTable(i,2));
                    }

                    parentDuct->pAcinus->VA  *= modTable(i,2);
                    parentDuct->pAcinus->VA0 *= modTable(i,2);

                    // Modify acinus stretch width
                    parentDuct->pAcinus->phi_dV = modTable(i,3);

                    // Modify acinus resistance
                    parentDuct->pAcinus->Racin *= modTable(i,4);
                }
            }
        }
        return;
    }

    // Call next generation
    applyModifications(parentDuct->pMajorDaughter);
    applyModifications(parentDuct->pMinorDaughter);

}


// Correct total volume after modification of ducts and acini
//**************************************************************/
void lung::correctTotalVolumeAfterModification(duct* parentDuct){

    // Declarations
    bool wasModified;

    // In ending duct apply correction
    if (parentDuct->reachedEnd(dLimit)){

        // Check whether the acinus volume has been modifed (then must not be corrected)
        if (nbrLinesModTab > 0){

            for (int i = 0; i<nbrLinesModTab; i++){
                if (parentDuct->pAcinus->getAcinusAbsIndex() == modTable(i,0)){


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

        // Correct the remaining acini
        if (!wasModified){
            parentDuct->pAcinus->VA  += TAVRed/(nbrAcini - nbrAciniRed);
            parentDuct->pAcinus->VA0 += TAVRed/(nbrAcini - nbrAciniRed);
        }
        TAV_aft_mod += parentDuct->pAcinus->VA;
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

        // transfer to inlet flow / velocity in acinus
        parentDuct->pAcinus->Q = parentDuct->Q;
        parentDuct->pAcinus->u = parentDuct->u;

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


// Update flow properties in acinus
//**************************************************************/
void lung::updateAcinus(duct* parentDuct){

    // Declarations
    double Flow, volumeRatio;

    // Update acinus volume when end duct is reached
    if (parentDuct->reachedEnd(dLimit)){

        // Set and get ratio of updated volume to former volume
        Flow = parentDuct->Q;
        parentDuct->pAcinus->setVolumeRatio(dt/Nt, Flow);
        volumeRatio = parentDuct->pAcinus->getVolumeRatio();

        // Update inlet flow and volume of trumpet acinus
        parentDuct->pAcinus->updateTrumpetAcinus(washout, dt/Nt, Flow);
        return;
    }

    // Call next generation
    updateAcinus(parentDuct->pMajorDaughter);
    updateAcinus(parentDuct->pMinorDaughter);
}


// Time integration of concentration in ducts and acini
//**************************************************************/
void lung::updateConcentrationInDucts(duct* parentDuct){

    // Declarations
    bool iT, rE;
    double cCondAW0, cCondAW1, c1Acinus0, c2Acinus0, c1Acinus1, c2Acinus1, dcdxCondAW0, dcdxCondAW1, hCondAW0, hCondAW1;

    // Pick arguments for different cases
    if (parentDuct->gen==0){
        iT = true;
        rE = false;
        c1Acinus0 = c2Acinus0 = c1Acinus1 = c2Acinus1 = -1.; // not used
    }
    else if (parentDuct->reachedEnd(dLimit)){
        iT = false;
        rE = true;
        if (washout==0){
            c1Acinus0 = parentDuct->pAcinus->pSpecies0->c1A;
            c2Acinus0 = parentDuct->pAcinus->pSpecies0->c2A;
        }
        if (washout==1){
            c1Acinus0 = parentDuct->pAcinus->pSpecies0->c1A;
            c2Acinus0 = parentDuct->pAcinus->pSpecies0->c2A;
            c1Acinus1 = parentDuct->pAcinus->pSpecies1->c1A;
            c2Acinus1 = parentDuct->pAcinus->pSpecies1->c2A;
        }
    }
    else{
        iT = false;
        rE = false;
        c1Acinus0 = c2Acinus0 = c1Acinus1 = c2Acinus1 = -1.; // not used
    }

    // Update concentration in current duct
    if (washout==0){
        parentDuct->pSpecies0->updateConcentrationInDuct(iT, rE, dt/Nt, cin, c1Acinus0, c2Acinus0);
    }
    if (washout==1){
        parentDuct->pSpecies0->updateConcentrationInDuct(iT, rE, dt/Nt, cin, c1Acinus0, c2Acinus0);
        parentDuct->pSpecies1->updateConcentrationInDuct(iT, rE, dt/Nt, cin, c1Acinus1, c2Acinus1);
    }

    // Anchor
    if (parentDuct->reachedEnd(dLimit)){

        // Define values to be passed to acinus
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

        // Update concentration in trumpet acinus
        if (washout==0){
            parentDuct->pAcinus->pSpecies0->updateConcentrationInTrumpetAcinus(dt/Nt, cCondAW0, dcdxCondAW0, hCondAW0);
        }
        if (washout==1){
            parentDuct->pAcinus->pSpecies0->updateConcentrationInTrumpetAcinus(dt/Nt, cCondAW0, dcdxCondAW0, hCondAW0);
            parentDuct->pAcinus->pSpecies1->updateConcentrationInTrumpetAcinus(dt/Nt, cCondAW1, dcdxCondAW1, hCondAW1);
        }
        return;
    }

    // Call next generation
    updateConcentrationInDucts(parentDuct->pMajorDaughter);
    updateConcentrationInDucts(parentDuct->pMinorDaughter);

}


// Collect data in trumpet acinus
//**************************************************************/
void lung::collectAllOutputDataInTrumpetAcinus(duct* parentDuct){

    // Declarations
    int N, fni, fzi, fcri, fcni, k2, k3;
    double x, y, z, dx, dy, dz, dp;
    double dist, lA, L, a;

    N = parentDuct->pAcinus->pSpecies0->N;
    fni = firstNodeIndA;
    fzi = firstZoneIndA;
    fcri = firstCoordIndA;
    fcni = firstConnIndA;

    // Node coordinates
    for (int k=0; k<=N; k++){
        k3 = 3*k;

        lA = parentDuct->pAcinus->lA;
        L  = parentDuct->pAcinus->L;

        dist = k*lA/double(N) / 5.;
        x = parentDuct->x[1] + cos(parentDuct->phi_z)*cos(parentDuct->phi)*dist;
        y = parentDuct->y[1] + cos(parentDuct->phi_z)*sin(parentDuct->phi)*dist;
        z = parentDuct->z[1] + sin(parentDuct->phi_z)*dist;


        pCoordinatesA[fcri+k3]   = x;
        pCoordinatesA[fcri+k3+1] = y;
        pCoordinatesA[fcri+k3+2] = z;
    }

    // Zone types
    for (int k=0; k<N; k++){
        pZoneTypesA[fzi+k] = 3; // is VISIT_LINE (see visit_writer.h)
    }

    // Connectivity
    for (int k=0; k<N; k++){
        k2 = 2*k;

        pConnectivityA[fcni+k2] = fni + k;
        pConnectivityA[fcni+k2+1] = fni + k + 1;
    }


    // Nodal data
    for (int k=0; k<=N; k++){

        // Absolut index
        pNodalAbsIndA[fni+k] = parentDuct->pAcinus->indAbsAcin;

        // Concentration values
        if (washout==0){
            pNodalConcentrationIA[fni+k]  = parentDuct->pAcinus->pSpecies0->C(k);
            pNodalConcentrationIIA[fni+k] = -1;
        }
        if (washout==1){
            pNodalConcentrationIA[fni+k]  = parentDuct->pAcinus->pSpecies0->C(k);
            pNodalConcentrationIIA[fni+k] = parentDuct->pAcinus->pSpecies1->C(k);
        }

        if (parentDuct->pAcinus->pSpecies0->C(k) > 1.01){
          cout << "ERROR: Concentration in trumpet acinus higher than 1: C(k) = " << parentDuct->pAcinus->pSpecies0->C(k) << endl;
        }

        // Radius
        pNodalRadiusA[fni+k] = 0.1*sqrt(4. * parentDuct->pAcinus->pSpecies0->ATr(k)/pi)/2.;

        // Velocity
        pNodalVelocityA[fni+k] = parentDuct->pAcinus->pSpecies0->UphTr(k);

        // Acinus volume
        pNodalAcinusVolumeA[fni+k] = parentDuct->pAcinus->VA - parentDuct->pAcinus->VA0;

        // Pleural pressure
        pNodalPleuralPressureA[fni+k] = ppl;
    }


    // Increment indices;
    firstNodeIndA += (N + 1);
    firstZoneIndA +=  N;
    firstCoordIndA += 3*(N + 1);
    firstConnIndA += 2*N;

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
        pNodalAbsIndD[fni+k] = parentDuct->indAbs;

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

        // Radius
        pNodalRadiusD[fni+k] = parentDuct->d/2.;


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

        // Collect in trumpet acinus
        collectAllOutputDataInTrumpetAcinus(parentDuct);

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
    firstNodeIndA = firstZoneIndA = firstCoordIndA = firstConnIndA = 0;
    collectAllOutputData(parentDuct);

    // Write duct data
    pAllOutputVariablesD[0] = pNodalAbsIndD;
    pAllOutputVariablesD[1] = pZonalVelocityD;
    pAllOutputVariablesD[2] = pNodalPressureD;
    pAllOutputVariablesD[3] = pNodalConcentrationID;
    pAllOutputVariablesD[4] = pNodalConcentrationIID;
    pAllOutputVariablesD[5] = pNodalRadiusD;
    pAllOutputVariablesD[6] = isModifiedD;


    char filenameD[100];
    sprintf(filenameD, "data/duct/TotalLungOutputDuct%05d.vtk", frame);

    write_unstructured_mesh(filenameD, 0, NntotD, pCoordinatesD, NztotD, pZoneTypesD, pConnectivityD, NvarD, pDimVarD, pCenteringD, pVarnamesD, pAllOutputVariablesD);

    // Write acinus data
    pAllOutputVariablesA[0] = pNodalAbsIndA;
    pAllOutputVariablesA[1] = pNodalVelocityA;
    pAllOutputVariablesA[2] = pNodalConcentrationIA;
    pAllOutputVariablesA[3] = pNodalConcentrationIIA;
    pAllOutputVariablesA[4] = pNodalRadiusA;
    pAllOutputVariablesA[5] = pNodalAcinusVolumeA;
    pAllOutputVariablesA[6] = pNodalPleuralPressureA;
    pAllOutputVariablesA[7] = isModifiedA;

    char filenameA[100];
    sprintf(filenameA, "data/acinus/TotalLungOutputAcinus%05d.vtk", frame);

    write_unstructured_mesh(filenameA, 0, NntotA, pCoordinatesA, NztotA, pZoneTypesA, pConnectivityA, NvarA, pDimVarA, pCenteringA, pVarnamesA, pAllOutputVariablesA);
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

#include "duct.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

using namespace std;


// Constructor
//**************************************************************/
duct::duct(controlProperties *contProp_, systemProperties *sysProp_){

    // Write class members
    //**********************************************************/
    contProp            = contProp_;
    sysProp             = sysProp_;

    // Initial variables
    endCAW = false;
    gen = ind = -1;
    indAbs = 0;
    d = l = lc = VD = 0;

    scal = 0;

    T = K1 = K2 = R = 0.;
    p = 101325;
    u = 0;
    Q = 0;

    fmod = 1.;

    bilinIntA = VectorXd::Zero(4);
    bilinIntC = VectorXd::Zero(4);
    bilinIntM = MatrixXd::Zero(4,4);

    pParent = pMinorDaughter = pMajorDaughter = 0;
    pLobule = 0;
    pSpecies0 = pSpecies1 = 0;


}

// Destructor
//**************************************************************/
duct::~duct(){
    // delete pMinorDaughter;
    // delete pMajorDaughter;

    // delete pLobule;
}


// Write Data (currently unused)
//**************************************************************/
void duct::writeData(FILE* aFile){

    fprintf(aFile, "%i %i %i %f %f %f %f %f %f %f %i %f\n",
            indAbs, gen, ind, d, l, x[0], x[1], y[0], y[1], phi, endCAW, scal);

}


// Check whether limit is reached
//**************************************************************/
bool duct::reachedEnd(double dLimit){

    if (d <= dLimit){
        endCAW = true;
        pMinorDaughter = pMajorDaughter = 0;
        return true;
    }
    else{
        return false;
    }
}


// Grow duct with daughters
//**************************************************************/
void duct::grow(int indA, double* kappaA, double* Phi){


    if (indA > 1 || indA < 0){
        cout << "ERROR: index for daughter duct has to be 0 or 1." << endl;
    }

    // Add new duct for daughter
    duct* daughter = new duct(contProp, sysProp);

    // Define properties
    daughter->gen = gen+1;
    daughter->ind = 2*ind + indA + 1;

    // Take data according to Weibel until Segmented Bronchus
    switch (daughter->gen){
        case 1:
            daughter->d = 0.68*d;
            daughter->l = 0.40*l*0.5; // compensate for larynx and mounth and apperatus
            break;

        case 2:
            daughter->d = 0.68*d;
            daughter->l = 0.40*l;
            break;
        case 3:
            daughter->d = 0.67*d;
            daughter->l = 0.40*l;
            break;

        case 4:
            daughter->d = 0.80*d;
            daughter->l = 1.67*l;
            break;
    }

    if (daughter->gen > 4){
        daughter->d   = kappaA[indA]*d;
        daughter->l   = kappaA[indA]*l;
    }

    // Cumulative length
    daughter->lc = lc + kappaA[indA]*l;

    // Coordinates for visualization
    daughter->phi  = Phi[indA]+phi;
    daughter->phi_z  = Phi[indA+2]+phi_z;
    daughter->x[0] = x[1];
    daughter->x[1] = x[1] + daughter->l*cos(daughter->phi_z)*cos(daughter->phi);
    daughter->y[0] = y[1];
    daughter->y[1] = y[1] + daughter->l*cos(daughter->phi_z)*sin(daughter->phi);
    daughter->z[0] = z[1];
    daughter->z[1] = z[1] + daughter->l*sin(daughter->phi_z);

    // Relate daughter to this parent
    daughter->pParent = this;

    switch (indA) {
        case 0:
            pMajorDaughter = daughter;
            break;

        case 1:
            pMinorDaughter = daughter;
            break;

        default:
            break;
    }

}


// Add species in duct
//**************************************************************/
void duct::addSpecies(int sp, int wo, double gridsize){

    if (wo==0) { // MBW with one species definied by 'sp'
        gas* species0 = new gas(contProp, sysProp, true, false, sp, wo, gridsize, l);

        pSpecies0 = species0;

        // Diffusion cross section
        pSpecies0->setCrossSection(d);
    }

    else if (wo==1) { // DTG-SBW
        gas* species0 = new gas(contProp, sysProp, true, false, 2, wo, gridsize, l); // helium (He)
        gas* species1 = new gas(contProp, sysProp, true, false, 3, wo, gridsize, l); // sulfur hexafluorid (SF6)

        pSpecies0 = species0;
        pSpecies1 = species1;

        // Diffusion cross section
        pSpecies0->setCrossSection(d);
        pSpecies1->setCrossSection(d);
    }
    else {
        cout << "ERROR: washout index is not known" << endl;
        return;
    }
}


// Relate gas classes from connecting ducts to each other
void duct::relateSpecies(int wo){

    if (wo==0) { // MBW with one species definied by 'sp'
        if (gen>0){
            pSpecies0->pParent = pParent->pSpecies0;
        }

        if (!endCAW){
            pSpecies0->pMajorDaughter = pMajorDaughter->pSpecies0;
            pSpecies0->pMinorDaughter = pMinorDaughter->pSpecies0;

            pSpecies0->pMajorDaughter->pSister = pMinorDaughter->pSpecies0;
            pSpecies0->pMinorDaughter->pSister = pMajorDaughter->pSpecies0;
        }
    }

    else if (wo==1) { // DTG-SBW
        if (gen>0){
            pSpecies0->pParent = pParent->pSpecies0;
            pSpecies1->pParent = pParent->pSpecies1;
        }

        if (!endCAW){
            pSpecies0->pMajorDaughter = pMajorDaughter->pSpecies0;
            pSpecies0->pMinorDaughter = pMinorDaughter->pSpecies0;
            pSpecies0->pMajorDaughter->pSister = pMinorDaughter->pSpecies0;
            pSpecies0->pMinorDaughter->pSister = pMajorDaughter->pSpecies0;

            pSpecies1->pMajorDaughter = pMajorDaughter->pSpecies1;
            pSpecies1->pMinorDaughter = pMinorDaughter->pSpecies1;
            pSpecies1->pMajorDaughter->pSister = pMinorDaughter->pSpecies1;
            pSpecies1->pMinorDaughter->pSister = pMajorDaughter->pSpecies1;
        }
    }
    else {
        cout << "ERROR: washout index is not known" << endl;
        return;
    }
}


// Set duct volume
//**************************************************************/
void duct::setDuctVolume(){

    double pi = 4.*atan(1.);

    VD = (d*d)/4.*pi*l;

}


// Return duct volume
//**************************************************************/
double duct::getDuctVolume(){

    if (VD == 0){
        cout << "NON-PHYSIOLOGICAL VALUE: duct volume is zero" << endl;
        return 0;
    }

    return VD;
}


// Connect lobule to duct
//**************************************************************/
void duct::connectLobule(bool scalingLbL, lobule* pLbT){

    // Create root lobule compartment and size with diameter
    lobule* rootLobule = new lobule(contProp, sysProp);

    // Initialize from template
    // Set lobule volume
    if (scalingLbL){
        rootLobule->VLb0 = pLbT->VLb0*pow(lc,1);
        rootLobule->VLb  = pLbT->VLb0*pow(lc,1);
    }
    else{
        rootLobule->VLb0 = pLbT->VLb0;
        rootLobule->VLb  = pLbT->VLb0;
    }

    // Set diameter of lobule
    rootLobule->d   = d;

    // Relate with feeding duct
    pLobule = rootLobule;
}


// Set index of lobule
//**************************************************************/
void duct::setDuctAbsIndex(int nbrDucts){

    if (nbrDucts < 0){
        cout << "ERROR: duct absolute index is smaller than zero" << endl;
    }

    indAbs = nbrDucts;
}


// Return index of lobule
//**************************************************************/
int duct::getDuctAbsIndex(){

    if (indAbs < 0){
        cout << "NON-PHYSIOLOGICAL VALUE: duct absolute index is smaller than zero" << endl;
        return 0;
    }

    return indAbs;
}


// Set transmissibility in current duct
//**************************************************************/
void duct::setTransmissibility(double mu, double TF){

    double pi = 4.*atan(1.);

    if (d <= 0 || l <= 0 || mu <= 0){
        cout << "ERROR: duct geometry (d, l) or dynamic viscosity not given" << endl;
        return;
    }

    // Compute pulsatile transmissibility
    T = fmod*TF*(pi*pow(d/2.,4))/(8.*mu*l);

    R = 1./T;
}

// Bilinear interpolation from pre-computed and read transmissibility factor table
//**************************************************************/
double duct::bilinearInterpolation(double evalTB, double evalD, int nbrLinesTransFact, MatrixXd& transFact, MatrixXd& transFactDomainTB, MatrixXd& transFactDomainD){

    // declarations
    int ii, jj;
    double TBl, TBu, dl, du;
    double TF;

    // find neighboring points
    jj = 0;
    for (int j=0; j<nbrLinesTransFact-1; j++){
        if ((transFactDomainTB(0,j) <= evalTB) && (evalTB < transFactDomainTB(0, j+1))){
            jj  = j;
            TBl = transFactDomainTB(0,j);
            TBu = transFactDomainTB(0,j+1);
        }
    }

    ii = 0;
    for (int i=0; i<nbrLinesTransFact-1; i++){
        if ((transFactDomainD(i,0) <= evalD) && (evalD < transFactDomainD(i+1, 0))){
            ii  = i;
            dl = transFactDomainD(i,0);
            du = transFactDomainD(i+1,0);
        }
    }

    // construct linear system
    bilinIntC(0) = transFact(ii,jj);
    bilinIntC(1) = transFact(ii,jj+1);
    bilinIntC(2) = transFact(ii+1,jj);
    bilinIntC(3) = transFact(ii+1,jj+1);

    bilinIntM(0,0) = 1; bilinIntM(0,1) = TBl; bilinIntM(0,2) = dl; bilinIntM(0,3) = TBl*dl;
    bilinIntM(1,0) = 1; bilinIntM(1,1) = TBu; bilinIntM(1,2) = dl; bilinIntM(1,3) = TBu*dl;
    bilinIntM(2,0) = 1; bilinIntM(2,1) = TBl; bilinIntM(2,2) = du; bilinIntM(2,3) = TBl*du;
    bilinIntM(3,0) = 1; bilinIntM(3,1) = TBl; bilinIntM(3,2) = du; bilinIntM(3,3) = TBu*du;

    // solve for coefficients
    bilinIntA = bilinIntM.householderQr().solve(bilinIntC);

    // evaluate function
    TF = bilinIntA(0) + bilinIntA(1)*evalTB + bilinIntA(2)*evalD + bilinIntA(3)*evalTB*evalD;

    return TF;

}


// Get transmissibility of duct
//**************************************************************/
double duct::getTransmissibility(){

    if (T <= 0){
        cout << "MISSING VALUE: duct transmissibility not given" << endl;
        return 0;
    }

    return T;
}


// Set K1 for transmissibility in end ducts
//**************************************************************/
void duct::setK1(double dt){

    // Declarations
    double E, Racin;
    double VLb, VLb0, Vstar;
    double dP, dV, gamma, beta, theta;

    // Get material law coefficients and Resistance of lobule
    E = pLobule->E;
    Racin = pLobule->Racin;
    gamma = pLobule->gamma;

    // Get lobule stretch
  	dP = pLobule->dP;
  	dV = pLobule->dV;

    if (E <= 0){
        cout << "MISSING VALUE: elasticity of lobule not given" << endl;
    }

    if (Racin <= 0){
        cout << "MISSING VALUE: resistance of lobule not given" << endl;
    }
    if (dP <= 0){
        cout << "MISSING VALUE: pressure stretch not given" << endl;
    }

    if (dV <= 0){
        cout << "MISSING VALUE: volume stretch not given" << endl;
    }

    // get volumes
    VLb  = pLobule->VLb;
    VLb0 = pLobule->VLb0;

    // non-linear constitutive law
    beta = dP/(exp(gamma*(VLb0 + dV)) - exp(gamma*VLb0));

    theta = 0.8;

    Vstar = VLb + dt*Q;
    K1 = (Racin*T + theta*dt*T*beta*gamma*exp(gamma*Vstar));
}


// Set K2 for rhs
//**************************************************************/
void duct::setK2(int bc, int nbrEndDucts, double TV, double dt, double varinp){

    // inlet flow b. c.       (bc = 0) => varinp = ppl    (current pleural pressure, computed in preceeding step)
    // pleural pressure b. c. (bc = 1) => varinp = dppldt (time derivative of pleural pressure)

    // declarations
    double E, Racin;
    double VLb, VLb0, Vstar;
    double dP, dV, gamma, beta, theta;
    double pm1, pPm1;

    // get material law coefficients and Resistance of lobule
    E = pLobule->E;
    Racin = pLobule->Racin;
    gamma = pLobule->gamma;

    // get lobule stretch
  	dP = pLobule->dP;
  	dV = pLobule->dV;

    if (E <= 0){
        cout << "MISSING VALUE: elasticity of lobule not given" << endl;
    }

    if (Racin <= 0){
        cout << "MISSING VALUE: resistance of lobule not given" << endl;
    }
    if (dP <= 0){
        cout << "MISSING VALUE: pressure stretch not given" << endl;
    }

    if (dV <= 0){
        cout << "MISSING VALUE: volume stretch not given" << endl;
    }

    // get volumes
    VLb  = pLobule->VLb;
    VLb0 = pLobule->VLb0;

    // old pressure values are current pressure values
    pPm1 = pParent->p;
    pm1  = p;

    // non linear constitutive law
    beta = dP/(exp(gamma*(VLb0 + dV)) - exp(gamma*VLb0));

    theta = 0.8;

    Vstar = VLb + dt*Q;
    if (bc==0){
        K2 = (1. + Racin*T)*pm1 - Racin*T*pPm1    - varinp + (1.-theta)*dt*T*gamma*beta*exp(gamma*Vstar)*(pPm1 - pm1);
    }
    else{
        K2 = (1. + Racin*T)*pm1 - Racin*T*pPm1 - dt*varinp + (1.-theta)*dt*T*gamma*beta*exp(gamma*Vstar)*(pPm1 - pm1);
    }
}


// Return K1
//**************************************************************/
double duct::getK1(){
    return K1;
}


// Return K2
//**************************************************************/
double duct::getK2(){
    return K2;
}


// Write coefficients into matrix
//**************************************************************/
void duct::writeCoeff(int bc, int nbrDucts, MatrixXd& coeffMat){

    // Declarations
    int indAbsm1, indAbsp10, indAbsp11;
    double Tp10, Tp11;

    // For trachea
    if (indAbs == 0){

        // Get transmissibilty
        Tp10 = pMajorDaughter->getTransmissibility();
        Tp11 = pMinorDaughter->getTransmissibility();

        // Get absolute indices
        indAbsp10 = pMajorDaughter->getDuctAbsIndex();
        indAbsp11 = pMinorDaughter->getDuctAbsIndex();

        // Write entries in coefficient matrix
        if (bc==0){ // inlet flow boundary conditions => solve for pleural pressure
            coeffMat(indAbs+1, indAbs)    += -T - Tp10 - Tp11;
            coeffMat(indAbs+1, indAbsp10) +=  Tp10;
            coeffMat(indAbs+1, indAbsp11) +=  Tp11;

            // ADDITIONAL EQUATION from inlet flow condition
            coeffMat(indAbs, indAbs)  += -T;
        }

        else{ // pleural pressure boundary conditions
            coeffMat(indAbs, indAbs)    += -T - Tp10 - Tp11;
            coeffMat(indAbs, indAbsp10) +=  Tp10;
            coeffMat(indAbs, indAbsp11) +=  Tp11;
        }

    }

    // For end ducts
    else if (endCAW){

        // get absolute indices
        indAbsm1 = pParent->getDuctAbsIndex();

        // write entries in coefficient matrix
        if (bc==0){ // inlet flow boundary conditions => solve for pleural pressure
            coeffMat(indAbs+1, indAbsm1) -= K1;
            coeffMat(indAbs+1, indAbs)   += (1+K1);
            // entry for pleural pressure
            coeffMat(indAbs+1, nbrDucts) -=  1;
        }

        else{ // pleural pressure boundary conditions
            coeffMat(indAbs, indAbsm1) -= K1;
            coeffMat(indAbs, indAbs)   += (1+K1);
        }
    }

    // For rest of airway
    else {

        // Get transmissibility
        Tp10 = pMajorDaughter->getTransmissibility();
        Tp11 = pMinorDaughter->getTransmissibility();

        // Get absolute indices
        indAbsm1 = pParent->getDuctAbsIndex();
        indAbsp10 = pMajorDaughter->getDuctAbsIndex();
        indAbsp11 = pMinorDaughter->getDuctAbsIndex();

        // Wire entries in coefficient matrix
        if (bc==0){ // inlet flow boundary conditions => solve for pleural pressure
            coeffMat(indAbs+1, indAbsm1)  +=  T;
            coeffMat(indAbs+1, indAbs)    += -T - Tp10 - Tp11;
            coeffMat(indAbs+1, indAbsp10) += Tp10;
            coeffMat(indAbs+1, indAbsp11) += Tp11;
        }

        else{ // pleural pressure boundary conditions
            coeffMat(indAbs, indAbsm1)  +=  T;
            coeffMat(indAbs, indAbs)    += -T - Tp10 - Tp11;
            coeffMat(indAbs, indAbsp10) += Tp10;
            coeffMat(indAbs, indAbsp11) += Tp11;
        }
    }
}


// Write rhs entry for each duct
//**************************************************************/
void duct::writeRHS(int bc, int nbrEndDucts, double TV, double dt, double patm, double varinp, double Qin, VectorXd& rhsVec){

    // Declarations
    double pP;

    //For trachea
    if (indAbs == 0){

        if (bc==0){ // inlet flow boundary conditions

            // Boundary condition: pressure at the mouth
            rhsVec(indAbs+1) = -T*patm;

            // ADDITIONAL EQUATION
            // Boundary condition: inlet flow at mouth
            rhsVec(indAbs) = Qin - T*patm;
        }

        else{ // pleural pressure boundary conditions

            // Boundary condition: pressure at the mouth
            rhsVec(indAbs) = -T*patm;
        }

    }

    // For end ducts
    else if (endCAW){
        // boundary condition: flux into lobule compartment
        setK2(bc, nbrEndDucts, TV, dt, varinp);

        if (bc==0){ // inlet flow boundary conditions
            rhsVec(indAbs+1) = getK2();
        }

        else{ // pleural pressure boundary conditions
            rhsVec(indAbs) = getK2();
        }
    }

    // For rest of airways
    else{
        return;
    }

}


// Calculate flow in duct
//**************************************************************/
void duct::setFlow(int wo, double patm){

    double pi = 4.*atan(1.);

    // Declarations
    double pm1;

    // For trachea
    if (indAbs == 0){
        pm1 = patm;
    }

    // Elsewhere
    else{
        pm1 = pParent->p;
    }

    // Flow
    Q = T*(pm1-p);

    // Mean velocity
    u = Q/((d*d)/4.*pi);

    // Copy to species objects
    if (wo == 0){
        pSpecies0->u = u;
        pSpecies0->Q = Q;
        pSpecies0->computeEffectiveDiffCoeff(d);
        pSpecies0->setDiffusionCoefficient(1); // set effective Diffusion coefficient after Taylor dispersion
    }
    if (wo == 1){
        pSpecies0->u = u;
        pSpecies1->u = u;
        pSpecies0->Q = Q;
        pSpecies1->Q = Q;
        pSpecies0->computeEffectiveDiffCoeff(d);
        pSpecies1->computeEffectiveDiffCoeff(d);
        pSpecies0->setDiffusionCoefficient(1); // set effective Diffusion coefficient after Taylor dispersion
        pSpecies1->setDiffusionCoefficient(1); // set effective Diffusion coefficient after Taylor dispersion
    }
}


// Estimate flow in current duct
//**************************************************************/
void duct::estimateFlow(int wo, double patm, double Qin){

    double pi = 4.*atan(1.);

    // Declarations
    double pm1;

    // For trachea
    if (indAbs == 0){
        pm1 = patm;
    }

    // Elsewhere
    else{
        pm1 = pParent->p;
    }

    // Flow
    Q = T*(pm1-p);

    // Mean velocity
    u_max = Q/((d*d)/4.*pi);
}

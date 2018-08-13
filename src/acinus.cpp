#include "acinus.h"
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

using namespace std;


// Constructor
//**************************************************************/
acinus::acinus(controlProperties *contProp_, systemProperties *sysProp_, transportProperties *transProp_){


    // Write class members
    //**********************************************************/
    contProp            = contProp_;
    sysProp             = sysProp_;
    transProp           = transProp_;

    // Read properties from input structs
    lAm_fac = sysProp->lAm_fac;
    bending = sysProp->bending;
    kappa = sysProp->kappa;
    kappa_hat = 2.*kappa*kappa;

    // Initialize parameter
    E = Racin = 0;

    indAbsAcin = -1;

    phi_dV = 1.;
    dP = dV = 0;

    gamma = 5000000;

    VA = VA0 = Vtilde = tACV = 0;

    d = 0;

    lt = lA = 0;

    f = 0;

    r = 1;

    R = p = Q = u = 0;

    Qold = 0;

    pSpecies0 = pSpecies1 = 0;

    //     pFeedDuct = 0;

}


// Destructor
//**************************************************************/
acinus::~acinus(){

    //     delete pO0;
    //     delete pO1;
    //
    //     delete pH0;
    //     delete pH1;
    //     delete pH2;
    //     delete pH3;

}


// Set acinus index
//**************************************************************/
void acinus::setAcinusAbsIndex(int nbrAcini){

    indAbsAcin = nbrAcini;
}


// Get acinus index
//**************************************************************/
int acinus::getAcinusAbsIndex(){

    if (indAbsAcin < 0){
        cout << "MISSING VALUE: acinus absolute index has not yet been set" << endl;
        return 0;
    }

    return indAbsAcin;
}


// Set length of trumpet
void acinus::setTrumpetLength(bool scalingTL, double scaleTL){

    // define trumpet lobule length either scaled based on (mean) cumulative duct length or unscaled (same length for all trumpets)
    if (scalingTL){
        lA = scaleTL*pow(lAm_fac*VA0, 1./3.);
    }
    else{
        lA = pow(lAm_fac*VA0, 1./3.);
    }

    // derive "terminal duct length" for trumpet shape constraints from trumpet lobule length
    lt = lA * (1. - kappa)/(0.98*kappa);
    //cout << "lt = " << lt << "; lA = " << lA << "; VA = " << VA << endl;
}


// Add gas species in trumpet
//**************************************************************/
void acinus::addSpeciesInTrumpet(int sp, int wo, double terminalDuctGridSize){

    if (wo==0) { // MBW with one species defined by 'sp'
        gas* species0 = new gas(contProp, sysProp, false, true, sp, wo, terminalDuctGridSize, lA);

        pSpecies0 = species0;

        // Inlet cross-section for trumpet lobule (same as cross-section of terminal duct)
        pSpecies0->setCrossSection(d);

        // Cross section for trumpet
        pSpecies0->setTrumpetTransportDomain(lt, VA0);

        // set flage for debug promts
        if (indAbsAcin == 10){
            pSpecies0->promptFlag = true;
        }
        else{
            pSpecies0->promptFlag = false;
        }

    }

    else if (wo==1) { // DTG-SBW
        gas* species0 = new gas(contProp, sysProp, false, true, 2, wo, terminalDuctGridSize, lA); // helium (He)
        gas* species1 = new gas(contProp, sysProp, false, true, 3, wo, terminalDuctGridSize, lA); // sulfur hexafluorid (SF6)

        pSpecies0 = species0;
        pSpecies1 = species1;

        // Cross-section for last duct
        pSpecies0->setCrossSection(d);
        pSpecies1->setCrossSection(d);

        // Cross section for trumpet
        pSpecies0->setTrumpetTransportDomain(lt, VA0);
        pSpecies1->setTrumpetTransportDomain(lt, VA0);

    }

    else {
        cout << "ERROR: washout index is not known" << endl;
        return;
    }
}


// Get volume of acinus
//**************************************************************/
double acinus::getAcinusVolume(){

    if (VA == 0){
        cout << "NON-PHYSIOLOGICAL VALUE: acinus volume is zero" << endl;
        return 0;
    }

    return VA;
}


// Calculate Stretch Distribution
//**************************************************************/
void acinus::calculateStretchDistribution(int nbrAcini){

    // Declarations
    double stretchWidth, lowerDistLim, upperDistLim;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    stretchWidth = 0.05;
    lowerDistLim = 0.5;
    upperDistLim = 3.0;

    default_random_engine generator(seed);

    //lognormal_distribution<double> distribution(0.0, 0.2);
    //normal_distribution<double> distribution(1.0, 0.1);
    //uniform_real_distribution<double> distribution(1.0-stretchWidth, 1.0+stretchWidth);

    if (nbrAcini <= 0){
        cout << "MISSING VALUE: total number of acini not given" << endl;
    }

    if (indAbsAcin < 0){
        cout << "MISSING VALUE: acinus absolute index not yet defined" << endl;
    }

    // Distribution of volume difference parameter (generic, regular, with stretchWidth) ...
    //phi_dV = -2*stretchWidth/nbrAcini*indAbsAcin + 1 + stretchWidth;

    // ... or by a stochastic ditribution distribution as defined above
    //phi_dV = distribution(generator);

    // ... or just a constant normal value of '1' for all trumpets (and maybe further changes throught mod. file inputs)
    phi_dV = 1.0;

    if (phi_dV < lowerDistLim){
      phi_dV = lowerDistLim;
    }

    if (phi_dV > upperDistLim){
      phi_dV = upperDistLim;
    }

    //cout << phi_dV << endl;
}


// Calculate Stiffness of Acinus
//**************************************************************/
void acinus::calculateStiffnessParameters(int nbrAcini, double TV){

    // Declarations
    double dP_tilde, dV_tilde;

    if (nbrAcini <= 0){
        cout << "MISSING VALUE: total number of acini not given" << endl;
    }

    if (TV <= 0){
        cout << "MISSING VALUE: tidal volume not given" << endl;
    }

    // Set stretch range
    dP = 1500;
    dV = phi_dV*TV/nbrAcini;

    // Set stretch middle-point
    dP_tilde = 0.5*dP*(1. - bending);
    dV_tilde = 0.5*dV*(1. + bending);

    // Set elasticity for linear model
    E = dP/dV;

    // Set shape parameter gamma for non-linear model
    gamma = findGamma(1e03, dP, dV, dP_tilde, dV_tilde);
}



// Set volume ratio to previous time step
//**************************************************************/
void acinus::setVolumeRatio(double dt, double Flux){

    if (VA <= 0){
        cout << "NON-PHYSIOLOGICAL VALUE: acinus volume is not given." << endl;
        return;
    }

    // set ratio
    r = 1. + Flux*dt/VA;
}


// Get volume ratio to previous time step
//**************************************************************/
double acinus::getVolumeRatio(){

    if (r <= 0){
        cout << "MISSING VALUE: acinus volume is not yet computed." << endl;
    }

    return r;
}


// Update volume of trumpet acinus
//**************************************************************/
void acinus::updateTrumpetAcinus(int wo, double dt, double Flow){

    double pi = 4.*atan(1.);

    // Update inlet flow & velocity
    Q = Flow;
    u = 4.*Q/(pow(d,2)*pi);

    // Update volume
    // Euler method
    // VA += Q*dt;

    // Trapez rule
    VA += 0.5*(Q + Qold)*dt;
    Qold = Flow;

    // Change related properties trumpet acinus' species
    if (wo==0){
        pSpecies0->u = u;
        pSpecies0->Q = Q;
        pSpecies0->updateTrumpetProperties(VA);
        pSpecies0->computeEffectiveDiffCoeff(d);
        pSpecies0->setDiffusionCoefficient(1);
    }
    if (wo==1){
        pSpecies0->u = u;
        pSpecies0->Q = Q;
        pSpecies0->updateTrumpetProperties(VA);
        pSpecies0->computeEffectiveDiffCoeff(d);
        pSpecies0->setDiffusionCoefficient(1);
        pSpecies1->u = u;
        pSpecies1->Q = Q;
        pSpecies1->updateTrumpetProperties(VA);
        pSpecies1->computeEffectiveDiffCoeff(d);
        pSpecies1->setDiffusionCoefficient(1);
    }
}


// Find Gamma
//**************************************************************/
double acinus::findGamma(double g_0, double deltaP, double deltaV, double deltaP_tilde, double deltaV_tilde){

    // Settings
    int MAXIT = 10000;

    double RTOL = 1e-07;
    double ATOL = 1e-07;


    double g, gm1, F, dFdg;

    g = g_0;

    // Newton method
    for (int k=0; k<MAXIT; k++){

        if (k == MAXIT-1){
            cout << "NO-GOOD: maximum number of iterations reached in root-finding." << endl;
        }

        // Store
        gm1 = g;

        // Response of function
        F = deltaP_tilde - deltaP*(exp(g*deltaV_tilde) - 1.)/(exp(g*deltaV) - 1.);

        // Response of derivative of function
        dFdg = -deltaP*((deltaV_tilde*exp(deltaV_tilde*g))/(exp(deltaV*g) - 1.) - deltaV*exp(deltaV*g)*(exp(deltaV_tilde*g) - 1.)/pow((exp(deltaV*g) - 1.),2));

        // Iteration
        g -=  F/dFdg;

        if (abs(g - gm1) <= abs(g)*RTOL + ATOL){
            break;
        }
    }

    return g;

}

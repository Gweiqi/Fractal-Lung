#include "lobule.h"
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

using namespace std;


// Constructor
//**************************************************************/
lobule::lobule(controlProperties *contProp_, systemProperties *sysProp_){


    // Write class members
    //**********************************************************/
    contProp            = contProp_;
    sysProp             = sysProp_;

    // Read properties from input structs
    lLbm_fac = sysProp->lLbm_fac;
    bending = sysProp->bending;
    kappa = sysProp->kappa;
    kappa_hat = 2.*kappa*kappa;

    // Initialize parameter
    E = Racin = 0;

    indAbsLob = -1;

    phi_dV = 1.;
    dP = dV = 0;

    gamma = 5000000;

    VLb = VLb0 = Vtilde = tACV = 0;

    d = 0;

    lt = lLb = 0;

    f = 0;

    r = 1;

    R = p = Q = u = 0;

    Qold = 0;

    pSpecies0 = pSpecies1 = 0;

    //     pFeedDuct = 0;

}


// Destructor
//**************************************************************/
lobule::~lobule(){

    //     delete pO0;
    //     delete pO1;
    //
    //     delete pH0;
    //     delete pH1;
    //     delete pH2;
    //     delete pH3;

}


// Set lobule index
//**************************************************************/
void lobule::setLobuleAbsIndex(int nbrLobules){

    indAbsLob = nbrLobules;
}


// Get lobule index
//**************************************************************/
int lobule::getLobuleAbsIndex(){

    if (indAbsLob < 0){
        cout << "MISSING VALUE: lobule absolute index has not yet been set" << endl;
        return 0;
    }

    return indAbsLob;
}


// Set length of trumpet lobule
void lobule::setTrumpetLobuleLength(bool scalingLbL, double scaleLbL){

    // define trumpet lobule lobule length either scaled based on (mean) cumulative duct length or unscaled (same length for all trumpet lobules)
    if (scalingLbL){
        lLb = scaleLbL*pow(lLbm_fac*VLb0, 1./3.);
    }
    else{
        lLb = pow(lLbm_fac*VLb0, 1./3.);
    }

    // derive "terminal duct length" for trumpet lobule shape constraints from trumpet lobule length
    lt = lLb * (1. - kappa)/(0.98*kappa);
    //cout << "lt = " << lt << "; lLb = " << lLb << "; VLb = " << VLb << endl;
}


// Add gas species in trumpetlobule
//**************************************************************/
void lobule::addSpeciesInTrumpetLobule(int sp, int wo, double terminalDuctGridSize){

    if (wo==0) { // MBW with one species defined by 'sp'
        gas* species0 = new gas(contProp, sysProp, false, true, sp, wo, terminalDuctGridSize, lLb);

        pSpecies0 = species0;

        // Inlet cross-section for trumpet lobule (same as cross-section of terminal duct)
        pSpecies0->setCrossSection(d);

        // Cross section for trumpetlobule
        pSpecies0->setTrumpetLobuleTransportDomain(lt, VLb0);

        // set flage for debug promts
        if (indAbsLob == 10){
            pSpecies0->promptFlag = true;
        }
        else{
            pSpecies0->promptFlag = false;
        }

    }

    else if (wo==1) { // DTG-SBW
        gas* species0 = new gas(contProp, sysProp, false, true, 2, wo, terminalDuctGridSize, lLb); // helium (He)
        gas* species1 = new gas(contProp, sysProp, false, true, 3, wo, terminalDuctGridSize, lLb); // sulfur hexafluorid (SF6)

        pSpecies0 = species0;
        pSpecies1 = species1;

        // Cross-section for last duct
        pSpecies0->setCrossSection(d);
        pSpecies1->setCrossSection(d);

        // Cross section for trumpetlobule
        pSpecies0->setTrumpetLobuleTransportDomain(lt, VLb0);
        pSpecies1->setTrumpetLobuleTransportDomain(lt, VLb0);

    }

    else {
        cout << "ERROR: washout index is not known" << endl;
        return;
    }
}


// Get volume of lobule
//**************************************************************/
double lobule::getLobuleVolume(){

    if (VLb == 0){
        cout << "NON-PHYSIOLOGICAL VALUE: lobule volume is zero" << endl;
        return 0;
    }

    return VLb;
}


// Calculate Stretch Distribution
//**************************************************************/
void lobule::calculateStretchDistribution(int nbrLobules){

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

    if (nbrLobules <= 0){
        cout << "MISSING VALUE: total number of lobules not given" << endl;
    }

    if (indAbsLob < 0){
        cout << "MISSING VALUE: lobule absolute index not yet defined" << endl;
    }

    // Distribution of volume difference parameter (generic, regular, with stretchWidth) ...
    //phi_dV = -2*stretchWidth/nbrLobules*indAbsLob + 1 + stretchWidth;

    // ... or by a stochastic ditribution distribution as defined above
    //phi_dV = distribution(generator);

    // ... or just a constant normal value of '1' for all trumpetlobules (and maybe further changes throught mod. file inputs)
    phi_dV = 1.0;

    if (phi_dV < lowerDistLim){
      phi_dV = lowerDistLim;
    }

    if (phi_dV > upperDistLim){
      phi_dV = upperDistLim;
    }

    //cout << phi_dV << endl;
}


// Calculate Stiffness of lobule
//**************************************************************/
void lobule::calculateStiffnessParameters(int nbrLobules, double TV){

    // Declarations
    double dP_tilde, dV_tilde;

    if (nbrLobules <= 0){
        cout << "MISSING VALUE: total number of lobules not given" << endl;
    }

    if (TV <= 0){
        cout << "MISSING VALUE: tidal volume not given" << endl;
    }

    // Set stretch range
    dP = 1500;
    dV = phi_dV*TV/nbrLobules;

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
void lobule::setVolumeRatio(double dt, double Flux){

    if (VLb <= 0){
        cout << "NON-PHYSIOLOGICAL VALUE: lobule volume is not given." << endl;
        return;
    }

    // set ratio
    r = 1. + Flux*dt/VLb;
}


// Get volume ratio to previous time step
//**************************************************************/
double lobule::getVolumeRatio(){

    if (r <= 0){
        cout << "MISSING VALUE: lobule volume is not yet computed." << endl;
    }

    return r;
}


// Update volume of trumpet lobule
//**************************************************************/
void lobule::updateTrumpetLobule(int wo, double dt, double Flow){

    double pi = 4.*atan(1.);

    // Update inlet flow & velocity
    Q = Flow;
    u = 4.*Q/(pow(d,2)*pi);

    // Update volume
    // Euler method
    // VLb += Q*dt;

    // Trapez rule
    VLb += 0.5*(Q + Qold)*dt;
    Qold = Flow;

    // Change related properties in trumpet lobule species
    if (wo==0){
        pSpecies0->u = u;
        pSpecies0->Q = Q;
        pSpecies0->updateTrumpetLobuleProperties(VLb);
        pSpecies0->computeEffectiveDiffCoeff(d);
        pSpecies0->setDiffusionCoefficient(1);
    }
    if (wo==1){
        pSpecies0->u = u;
        pSpecies0->Q = Q;
        pSpecies0->updateTrumpetLobuleProperties(VLb);
        pSpecies0->computeEffectiveDiffCoeff(d);
        pSpecies0->setDiffusionCoefficient(1);
        pSpecies1->u = u;
        pSpecies1->Q = Q;
        pSpecies1->updateTrumpetLobuleProperties(VLb);
        pSpecies1->computeEffectiveDiffCoeff(d);
        pSpecies1->setDiffusionCoefficient(1);
    }
}


// Find Gamma
//**************************************************************/
double lobule::findGamma(double g_0, double deltaP, double deltaV, double deltaP_tilde, double deltaV_tilde){

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

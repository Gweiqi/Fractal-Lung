#ifndef ____acinus__
#define ____acinus__

#include "gas.h"
//#include "duct.h"
#include "dataStruct.h"

#include <cmath>
#include <iostream>

using namespace std;


/// Acinus
class acinus{

public:

    // Member Variables
    //**********************************************************/
    /// Pointer on control properties
    controlProperties *contProp;

    /// Pointer on system properties
    systemProperties *sysProp;

    /// Pointer on transport properties
    transportProperties *transProp;

    /// Absolut index of trumpet acinus
    int indAbsAcin;

    /// Stretch distribution
    double phi_dV;

    /// Strech range
    double dP, dV;

    /// Elasticity of acinus
    double E;

    /// Resistance of acinus
    double Racin;

    /// Non-linear stiffness parameter
    double gamma;

    /// Non-linear bending parameter (degree of non-linearity)
    double bending;

    /// Inlet diameter
    double d;

    /// Length of trumpet acinus
    double lA;

    /// Length of last conductive airway (lAm_fac from system properties, lt from scaling with true terminal duct)
    double lAm_fac, lt;

    /// Scaling parameter for length of trumpet acinus
    double lambda;

    /// Homothety rate of trumpet acinus
    double kappa_hat;

    /// Grow rate of trumpet acinus
    double kappa;

    /// flow factor in V-structure
    double f;

    /// limes for z (gen) -> inf of cummulative length
    double L;

    /// Volume of acinus, volume ratio (V-structure), relative volume
    double VA, r, Vtilde;

    /// Initial volume of trumpet acinus
    double VA0;

    /// Total alveolated cell volume
    double tACV;

    // variables for V-structure
    double R, p;

    /// Volume flow
    double Q;

    /// Volume flow of previous time step
    double Qold;

    /// Flow velocity
    double u;

    /// Stretch width for acini volume distribution
    double stretchWidth;

    /// Bool for modification
    bool isModified;

    /// Pointer on gas species 1
    gas* pSpecies0;

    /// Pointer on gas species 2
    gas* pSpecies1;

    /// Pointer on feeding duct
    //duct* pFeedDuct;

    // MEMBER FUNCTIONS
    //**********************************************************/
    /// Constructor of acinus
    acinus(controlProperties *contProp, systemProperties *sysProp, transportProperties *transProp);

    ~acinus();

    /// Set the absolute acinus index
    void setAcinusAbsIndex(int nbrAcini);

    /// Return the absolute acinus index
    int getAcinusAbsIndex();

    /// Add a class gas species in trumpet acinus
    void addSpeciesInTrumpet(int sp, int wo, double terminalDuctGridSize);

    /// Set the length of the trumpet
    void setTrumpetLength(bool scalingTL, double scaleTL);

    /// Calculate the acinus volume
    double getAcinusVolume();

    /// Calculate the volume distribution of the acini
	  void calculateStretchDistribution(int nbrAcini);

	  /// Set elasticity for the acinus
    void calculateStiffnessParameters(int nbrAcini, double TV);

    /// Set the volume ratio between input flux and acinus volume
    void setVolumeRatio(double dt, double Flux);

    /// Return volume ratio
    double getVolumeRatio();

    /// Update flow properties in the trumpet acinus
    void updateTrumpetAcinus(int wo, double dt, double Flow);

    /// Set shape parameter for non-linear model
    double findGamma(double g_0, double deltaP, double deltaV, double deltaP_tilde, double deltaV_tilde);

};


#endif /* defined(____acinus__) */

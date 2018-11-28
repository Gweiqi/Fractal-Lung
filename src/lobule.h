#ifndef ____lobule__
#define ____lobule__

#include "gas.h"
//#include "duct.h"
#include "dataStruct.h"

#include <cmath>
#include <iostream>

using namespace std;


/// Lobule
class lobule{

public:

    // Member Variables
    //**********************************************************/
    /// Pointer on control properties
    controlProperties *contProp;

    /// Pointer on system properties
    systemProperties *sysProp;

    /// Absolut index of trumpet lobule
    int indAbsLob;

    /// Stretch distribution
    double phi_dV;

    /// Strech range
    double dP, dV;

    /// Elasticity of lobule
    double E;

    /// Resistance of lobule
    double Racin;

    /// Non-linear stiffness parameter
    double gamma;

    /// Non-linear bending parameter (degree of non-linearity)
    double bending;

    /// Inlet diameter
    double d;

    /// Length of trumpet lobule
    double lLb;

    /// Length of last conductive airway (lLbm_fac from system properties, lt from scaling with true terminal duct)
    double lLbm_fac, lt;

    /// Scaling parameter for length of trumpet lobule
    double lambda;

    /// Homothety rate of trumpet lobule
    double kappa_hat;

    /// Grow rate of trumpet lobule
    double kappa;

    /// flow factor in V-structure
    double f;

    /// limes for z (gen) -> inf of cummulative length
    double L;

    /// Volume of lobule, volume ratio (V-structure), relative volume
    double VLb, r, Vtilde;

    /// Initial volume of trumpet lobule
    double VLb0;

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
    /// Constructor of lobule
    lobule(controlProperties *contProp, systemProperties *sysProp);

    ~lobule();

    /// Set the absolute lobule index
    void setLobuleAbsIndex(int nbrLobules);

    /// Return the absolute lobule index
    int getLobuleAbsIndex();

    /// Add a class gas species in trumpet lobule
    void addSpeciesInTrumpetLobule(int sp, int wo, double terminalDuctGridSize);

    /// Set the length of the trumpet
    void setTrumpetLobuleLength(bool scalingLbL, double scaleLbL);

    /// Calculate the lobule volume
    double getLobuleVolume();

    /// Calculate the volume distribution of the acini
	  void calculateStretchDistribution(int nbrLobules);

	  /// Set elasticity for the lobule
    void calculateStiffnessParameters(int nbrLobules, double TV);

    /// Set the volume ratio between input flux and lobule volume
    void setVolumeRatio(double dt, double Flux);

    /// Return volume ratio
    double getVolumeRatio();

    /// Update flow properties in the trumpet lobule
    void updateTrumpetLobule(int wo, double dt, double Flow);

    /// Set shape parameter for non-linear model
    double findGamma(double g_0, double deltaP, double deltaV, double deltaP_tilde, double deltaV_tilde);

};


#endif /* defined(____lobule__) */

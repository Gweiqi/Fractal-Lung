#ifndef ____duct__
#define ____duct__

#include "lobule.h"
#include "gas.h"
#include "dataStruct.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace Eigen;


/// Duct
class duct{

public:
    // Member Variables
    //**********************************************************/
    /// Pointer on control properties
    controlProperties *contProp;

    /// Pointer on system properties
    systemProperties *sysProp;

    /// Pointer on transport properties
    transportProperties *transProp;

    /// Bool for end of conducting airways
    bool endCAW;

    /// Generation
    int gen;

    /// Index of duct
    int ind;

    /// Absolute index of duct
    int indAbs;

    /// Diameter of duct
    double d;

    /// Length of duct
    double l;

    /// Cumulative length from trachea
    double lc;

    /// Volume of duct
    double VD;

    /// Scaling for lobules volume (only in terminal duct)
    double scal;

    // Visualization
    double phi, phi_z;
    double x[2];
    double y[2];
    double z[2];

    /// Pointer on parent duct
    duct* pParent;

    /// Pointer on minor daughter duct
    duct* pMinorDaughter;

    /// Pointer on major daughter duct
    duct* pMajorDaughter;

    /// Pointer on lobules compartment
    lobule* pLobule;

    /// Pointer on gas species 1
    gas* pSpecies0;

    /// Pointer on gas species 2
    gas* pSpecies1;

    /// Transmissibility of duct
    double T, R;

    /// Flux coefficients for end ducts
    double K1;

    /// Coeffient to build rhs
    double K2;

    /// Modification factor
    double fmod;

    /// Bool for modification
    bool isModified;

    /// Pressure in duct
    double p;

    /// Volume flux, mean velocity
    double Q;

    /// Velocity, maximum velocity
    double u, u_max;

    /// CFL number for diffusive term
    double CFL_Diff;

    /// CFL number for advective term
    double CFL_Adv;

    // data vectors for bilinear interpolation
    VectorXd bilinIntA, bilinIntC;
    MatrixXd bilinIntM;

    // MEMBER FUNCTIONS
    //**********************************************************/
    /// Constructor
    duct(controlProperties *contProp, systemProperties *sysProp, transportProperties *transProp);

    /// Destructor
    ~duct();

    /// Write data to file
    void writeData(FILE* aFile);

    /// Test for terminal bronchioles
    bool reachedEnd(double dLimit);

    /// Grow airway with one major/minor daughter
    void grow(int indA, double* kappaA, double* Phi);

    /// Add gas class to duct
    void addSpecies(int sp, int wo, double gridsize);

    /// Relate gass classes to parent, sister and daughters
    void relateSpecies(int wo);

    /// Compute duct volume
    void setDuctVolume();

    /// Return duct volume
    double getDuctVolume();

    /// Connect lobules class to end duct
    void connectLobule(bool scalingLbL, lobule* pLbT);

    /// Set the absolute index of duct
    void setDuctAbsIndex(int nbrDucts);

    /// Return absolute index of duct
    int getDuctAbsIndex();

    /// Calculate the transmissibility of current duct
    void setTransmissibility(double mu, double TF);

    // bilinear interpolation
    double bilinearInterpolation(double evalPointTB, double evalPointD, int nbrLinesTransFact, MatrixXd& transFact, MatrixXd& transFactDomainTB, MatrixXd& transFactDomainD);

    /// Return transmissibility
    double getTransmissibility();

    /// Set K1
    void setK1(double dt);

    /// Return K1
    double getK1();

    /// Set K2
    void setK2(int bc, int nbrEndDucts, double TV, double dt, double varinp);

    /// Get K2
    double getK2();

    /// Write transmissibility coefficients into matrix for pressure network
    void writeCoeff(int bc, int nbrDucts, MatrixXd& coeffMat);

    /// Write right-hand-side for pressure network
    void writeRHS(int bc, int nbrEndDucts, double TV, double dt, double patm, double varinp, double Qin, VectorXd& rhsVec);

    /// Set flow and velocity in duct after pressure network was solved
    void setFlow(int wo, double patm);

    /// Estimate flow to scale ducts with respect to the CFL condition
    void estimateFlow(int wo, double patm, double Qin);
};



#endif /* defined(____duct__) */

#ifndef ____gas__
#define ____gas__

#include <iostream>
#include <cmath>
#include <eigen_3_3_4/Dense>
#include <eigen_3_3_4/Sparse>
#include <eigen_3_3_4/IterativeLinearSolvers>

#include "dataStruct.h"

using namespace std;
using namespace Eigen;

// declare gas class
class gas{

public:

    /* VARIABLES
     -------------------*/
    // Pointer on control properties
    controlProperties *contProp;

    // Pointer on system properties
    systemProperties *sysProp;

    // Pointer on transport properties
    transportProperties *transProp;

    // flags
    bool chebychevcollocation;
    bool promptFlag;

    // species: 0_2 (0), CO_2 (1), He (2), SF_6 (3), N_2 (4)
    int species;

    // kind of washout: MBW (0) DTG SBW (1)
    int washout;

    // molar mass
    double mmass;

    // spatial discretization
    int N;
    double h;

    // trumpet acinus scaling and grow rate
    double kappa_hat, kappa;

    // modification factor for molecular diffusion coefficent
    double Diff_fac;

    // length and diameter of duct, trumpet
    double l, diam, lt, L;

    // advection
    double u, Q;

    // cross-sections ('St' at inlet of trumpet lobule)
    double S, St;

    // volume of trumpet acinus
    double VTr, VTrd;

    // diffusion coefficient
    double Dmol, Deff, Diff;
    MatrixXd DmolMat;

    // Peclet number
    double Pe;

    // mass flux
    double F_l, F_r;

    // interface surface
    double Sd_A, Sd_B, Sd_C;
    double fb_p, fb_maj, fb_min;

    // quantities at acinus inlet
    double c1A, c2A, cCAW, dcdxCAW;

    // time integration coefficient (Crank-Nicolson)
    double theta;

    // variables for chebychev collocation
    double zeta_i, zeta_j;
    double x_i, x_j;

    // streamwise coordinate
    double x;

    // trumpet geometry power law coefficients (p1, p2 for dynamic (time varing) trumpet lobule; p1d, p2d for advection (fixed size) trumpet)
    double p1, p2, p1d, p2d;

    // exponents used in power law trumpet
    double n1, n2;

    // intersection cross-section (constraint for power law)
    int z_star;
    double x_star, S_star;

    // pointers on relatives
    gas* pParent;
    gas* pMinorDaughter;
    gas* pMajorDaughter;
    gas* pSister;

    // species concentration (volumetric)
    VectorXd C;
    VectorXd Cm1;

    VectorXd rhs;

    // transport cross-section of trumpet acinus
    VectorXd STrd;

    // cross-section of trumpet lobule S = S(x,t)
    VectorXd ATr;

    // per-generation diameter of trumpet lobule
    VectorXd dTr;

    // physical velocity in trumpet acinus
    VectorXd UphTr;

    // finite difference and mapping operators operators
    MatrixXd I;

    MatrixXd D1_BW;
    MatrixXd D1_FW;
    MatrixXd D1_UW;
    MatrixXd D1_C, D2_C;

    MatrixXd D1_CC, D2_CC;

    // operators
    MatrixXd RO;
    MatrixXd LO;

    // effective advection and diffusion in trumpet acinus
    MatrixXd UTr;
    MatrixXd DIFFTr;

    // linear system for power law coefficients
    VectorXd powLawRhs;
    VectorXd powLawP;

    MatrixXd powLawM;

    // for quadratic interpolation
    VectorXd quadIntA; // coefficients of quadratic function (3)
    VectorXd quadIntC; // data points (3)
    VectorXd quadIntX; // control points (3)

    MatrixXd quadIntM; // matrix for linear problem involved in coefficient computation (3x3) (defined by control points)

    // iterative solver for concentration update
    Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverBiCGSTAB;

    /* MEMBER FUNCTIONS
     -------------------*/
    gas(controlProperties *contProp, systemProperties *sysProp_, bool condAW, bool trumpetacinus, int sp, int wo, double gridsize, double l);
    ~gas();

    void setCrossSection(double d);

    void setTrumpetTransportDomain(double terminalDuctLength, double VA0);

    void updateTrumpetProperties(double VA);

    void computeEffectiveDiffCoeff(double d);

    void setDiffusionCoefficient(int opt);

    void updateConcentrationInTrumpetAcinus(double dt,  double cLastCondAW,  double dcdxLastCondAW, double hLastCondAW);

    void updateConcentrationInDuct(bool inletTrachea, bool reachedEnd, double dt,  double cin, double c1Acinus, double c2Acinus);

    void updateOldConc();

    double quadraticInterpolation(double evalPoint, int opt1, int opt2);

};

#endif /* defined(____gas__) */

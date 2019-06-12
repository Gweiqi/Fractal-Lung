#ifndef ____lung__
#define ____lung__

#include "flDriver.h"
#include "duct.h"
#include "lobule.h"
#include "gas.h"
#include "IOdict.h"
#include "visit_writer.h"

#include <iostream>
#include <eigen/Dense>

using namespace std;
using namespace Eigen;


/// Lung
class lung{

public:
    // Variables
    //**********************************************************/
    /// Pointer on control properties
    controlProperties *contProp;

    /// Pointer on system properties
    systemProperties *sysProp;

    // flag for trumpet lobule scaling
    bool scalingLbL;

    // option for LPM boundary conditions (0: inlet flow rate. 1: pleural pressure)
    int bc;

    /// Maximum number of generation in the lobule
    int maxGenTrumpLob;

    /// Number of breath
    int nbrBreaths;

    /// Number of duct
    int nbrDucts;

    /// Number of end ducts
    int nbrEndDucts;

    /// Check number of end ducts
    int nbrEndDuctsCheck;

    /// Minimum generation
    int endDuctMinGen;

    /// Maximum Generation
    int endDuctMaxGen;

    /// Number of lobules
    int nbrLobules, nbrLobulesRed;

    /// Number of lines in modification table
    int nbrLinesModTab, nbrLinesTransFact;

    /// Kind of washout: MBW (0) DTG SBW (1)
    int washout;

    /// Species: 0_2 (0), CO_2 (1), He (2), SF_6 (3), N_2 (4)  - for MBW only
    int species;

    /// Number of grid points in smallest an longest duct
    int NgS, NgL;
    double umin;

    /// total number of grid points for data storage. Duct nodes; Duct zones; Lobule nodes; Lobule zones
    int NntotD, NztotD, NntotLb, NztotLb;

    // locatin of maximum required time step refinement
    int indAbs_Nt, gen_Nt;

    // pi
    double pi;

    /// total length of airways
    double Ltot;

    /// Inlet concentration
    double cin;

    /// Inlet volume flow
    double Qin;

    /// Outlet flow of conducting airways
    double Qouttot;

    /// Limit duct diameter
    double dLimit;

    /// Smallest duct length
    double lmin, dmin;

    /// Asymmetry parameters
    double r, eta;

    double kappaLb[2];
    double kappaH;
    double kappaV;

    /// Trumpet lobule grow rate
    double kappa, kappa_hat;

    /// Time integration parameters
    double dt, Tend, Nt;

    /// Dynamic viscosity
    double mu;

    /// Density
    double rho;

    /// Magnification factor for total lobular resistance
    double R_fac;

    /// Atmospheric pressure (pressure at the mouth)
    double patm;

    /// Diameter of trachea
    double dTrachea;

    /// Length of trachea
    double lTrachea;

    /// Pointer duct trachea
    duct* pTrachea;

    /// total airway volume
    double TAWV;

    /// Total lobule volume
    double TLbV, TLbV_aft_mod;

    /// Total residual volume (has to equal FRC)
    double TRLV;

    /// Total reduced lobule volume (due to modifications)
    double TLbVRed;

    /// Functional residual capacity; tidal volume; breath period
    double FRC;

    /// Tidal volume of current breath
    double TV;

    /// Mean tidal volume
    double TVm;

    /// Length of breath
    double TB;

    /// Mean length of breath
    double TBm;

    /// Matrix for tidal volume and breath length
    MatrixXd TBTV;

    /// Total lobule scaling
    double totalCompScaling;

    /// Lobule template root
    lobule* pLbTroot;

    /// Volume of lobule template
    double VLbT;

    /// Elasticity of lobule
    double E;

    // for visualization
    double Phi[4];

    /// Coefficient matrix for pressure network
    MatrixXd coeffMat;

    /// Pressure vector
    VectorXd pressureVec;

    /// Right-hand side vector
    VectorXd rhsVec;

    /// Inlet flow (input data)
    VectorXd inletFlow;
    int maxIndInletFlow;

    /// Pleura Pressure (input data)
    VectorXd pleuralPressure;

    /// Pulsatile transmissability factors for flow correction
    double d_min, d_max, TB_min, TB_max;
    MatrixXd transFact, transFactDomainD, transFactDomainTB;

    /// Table with modification values for ducts and lobules
    MatrixXd modTable;

    /// Pleural pressure
    double ppl, dppldt;

    /// simulation outputs (with visit_writer for VTK-i/o)
    IOdict writeConc;
    int frame;

    // Duct
    int* pZoneTypesD;
    int* pConnectivityD;

    float* pCoordinatesD;

    float* pNodalLbbsIndD;
    float* pZonalVelocityD;
    float* pNodalPressureD;
    float* pNodalConcentrationID;
    float* pNodalConcentrationIID;
    float* pNodalRadiusD;
    float*  isModifiedD;

    int NvarD;
    int* pDimVarD;
    int* pCenteringD;
    const char** pVarnamesD;

    float** pAllOutputVariablesD;

    int firstNodeIndD, firstZoneIndD, firstCoordIndD, firstConnIndD;

    // Trumpet lobule
    int* pZoneTypesLb;
    int* pConnectivityLb;

    float* pCoordinatesLb;

    float* pNodalLbbsIndLb;
    float* pNodalVelocityLb;
    float* pNodalConcentrationILb;
    float* pNodalConcentrationIILb;
    float* pNodalRadiusLb;
    float*  isModifiedLb;
    float* pNodalLbcinusVolumeLb;
    float* pNodalPleuralPressureLb;

    int NvarLb;
    int* pDimVarLb;
    int* pCenteringLb;
    const char** pVarnamesLb;

    float** pAllOutputVariablesLb;

    int firstNodeIndLb, firstZoneIndLb, firstCoordIndLb, firstConnIndLb;


    // Member functions
    //**********************************************************/
    /// Constructor
    lung(controlProperties *contProp, systemProperties *sysProp);

    /// Destructor
    ~lung();

    /// Deconstruct duct
    void deconstructLung(duct* parentDuct, duct* rubberD);

    /// Generate duct class of conducting airways
    void genMorphCondAirways(duct* parentDuct);

    /// Calculate conducting airway volume
    void airwayVolume(duct* parentDuct);

    /// Size lobules to match FRC
    void sizeLobuleTemplate();

    /// Attach lobules class to end ducts
    void genTrumpetLobuleOnAirways(duct* parentDuct);

    ///
    void lobuleVolumeNResistance(duct* parentDuct);

	  void calculateStiffnessParametersInTrumpetLobule(duct* parentDuct);

    /// Set gas class in ducts and lobules
    void setupGaseousSpecies(duct* parentDuct);

    /// Compute time step refinement in order to satisfy CFL condition in each element
    void computeTimeStepRefinement(duct* parentDuct, int opt);

    /// Compute total number of grid points
    void computeTotalNbrOfGridPoints(duct* parentDuct);

    /// Initialize matrix and vectors to solve pressure network
    void initializeMat();

    /// Initialize arrays for VTK output
    void initializeVTKArrays();

    /// Read TBTV data from file
    void readTBTV(breathFlow* allBreathFlow);

    /// Read inlet flow data from file
    void readInletFlow(breathFlow* allBreathFlow);

    /// Read pleural pressure data from file
    void readPleuralPressure(breathFlow* allBreathFlow);

    /// Read transition facts for pulsatile flow from file
    void readTransFact();

    /// Read modification table from file
    void readModifications();

    /// Apply modifications in ducts and lobules
    void applyModifications(duct* parentDuct);

    /// Correct volume after modification
    void correctTotalVolumeAfterModification(duct* parentDuct);

    /// Set TB and TV for current breath
    void setTBTV(int k);

    void airwaySetTK1(duct* parentDuct);

    /// Write coefficients in transmissibility matrix
    void writeCoeffMatrix(duct* parentDuct);

    // interpolate input data to sub-time step
    void interpInputData(int it, int n);

    /// Write right-hand-side for pressure network
    void writeRHSVector(duct* parentDuct, int opt);

    /// Distribute pressure after solving pressure network
    void distrPressure(duct* parentDuct);

    /// Compute flow with pressure values
    void computeFlow(duct* parentDuct);

    /// Estimate flow to scale spatial resolution (CFL condition)
    void estimateFlow(duct* parentDuct);

    /// Update condition in lobule
    void updateLobule(duct* parentDuct);

    /// Integrate advection diffusion equation in ducts and lobules
    void updateConcentrationInDucts(duct* parentDuct);

    /// Estimate resistance in lobules
    void computeTrumpetLobuleResistance(double d0, double l0, lobule* trumpetLobule);

    /// Collect output data in lobules for vtk output
    void collectAllOutputDataInTrumpetLobule(duct* parentDuct);

    /// Collect output data in ducts for vtk output
    void collectAllOutputData(duct* parentDuct);

    /// Write primary output data, concentration at inlet
    void writePrimaryLungData(int it, breathFlow* allBreathFlow, breathResults* allBreathResults, duct* parentDuct);

    /// Write full lung data, vtk files
    void writeFullLungData(duct* parentDuct);

    /// Write for debugging
    void writeOut(duct* parentDuct);

    /// Update old concentration in ducts
    void updateOldConc(duct* parentDuct);
};



#endif /* defined(____lung__) */

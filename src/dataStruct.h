#ifndef ____dataStruct__
#define ____dataStruct__

#include <fstream>
#include <vector>
#include <iostream>

using namespace std;


// Struct declarations
//**************************************************************/
///Struct containing data to initialize a lung object
struct init{
    double dt;              // Time stepping
    double Tend;            // Duration of simulation
    int washout;            // STG-MBW (0); MTG-SBW (1), defined in input file
    int species;            // species: 0_2 (0), CO_2 (1), He (2), SF_6 (3), N_2 (4)  - for MBW only
    int nbrBreaths;         // Number of breaths
    double dL;              // Limit for transitional bronchioles
    double FRC;             // Functional residual capacity FRC [in m^3]
    // double NgS, NgL;
};


/// Struct containing data for single breath
//**************************************************************/
struct breathFlow{
    int N;                  // Number of flow data points in this breath
    double* flowData;       // Sampled (and maybe re-sampled) flow values !!! in [m^3/s] !!!
    double* pplData;        // Generic pleural pressure data for option bc=1 (pleural pressure as boundary conditions)
    double TV;              // Tidal volume (average from expiration and inspiration) !!! in [m^3] !!!
    double TB;              // Breath periode (duration of one breath = inverse of breath rate) !!! in [s] !!!
    double fs;              // Sampling frequency of flow data !!! in [1/s] !!!
};


/// Struct containing conducting airway and lobules properties
//**************************************************************/
struct systemProperties{
    double cin;             // Inlet concentration
    int    bc ;             // Flag for boundary condition for LPM
    double p0;              // Ambient pressure
    double r;               // Bifurcation parameter
    double eta;             // Bifurcation parameter
    double kappa;           // Homothety in trumpet lobule
    int    maxGenLb;        // Number of generation in Lobules
    double lLbm_fac;        // mean trumpet lobule length factor
    double bending;         // bending of non-linear constitutive law for trumpet
    double n1,n2;           // Exponents used in power law for trumpet lobule
    double z_star;          // generation at which model cross-section intersects with exponential growth
    double R_fac;           // Magnification factor for total lobular resistance
    double Diff_fac;        // Modification factor for molecular diffusion coefficient
    int    scalingLbL;       // trumpet lobule scaling
};


/// Struct containing transport properties
//**************************************************************/
struct transportProperties{
    vector<double> Diff;    // Vector with diffusion coefficient for 1 carrier species [m^2/s]
    double rho;             // Density [kg/m^3]
    double mu;              // Dynamic viscosity [Pa*s]
};


/// Struct containing solver properties
//**************************************************************/
struct controlProperties{
    double dt;              // Time stepping
    double Tend;            // Duration of simulation
    int nbrBreaths;         // Number of breaths
    double generalCN;       // General Crank Nicolson coefficient
    double CFL;             // CFL number
    int NxDuctMin;          // Minimum number of grid points in ducts
    int NxDuctMax;          // Maximum number of grid points in ducts
    int NxLob;             // Number of grid points in acini
    int fOutFull;           // Write frequency full output
    bool writeFull;         // Write frequency full output (1 - yes, 0 - no)
    bool writeRedGrid;      // Full data output on reduced grid (1 - yes, 0 - no)
};


/// Struct containing solution data
//**********************************************************/
struct breathResults{ // contains results for one breath
    int N;                  // Number of data points in this breath

    // time stamp array
    double* time;

    // Computed species concentrations (normalized [0,1]);
    // DTG-SBW returns two species, therefore two arrays are given. In case of MBW simulation, the 2nd channel is filled with '-1' values.
    double* speciesIData; // Helium He
    double* speciesIIData; // Sulfur-Hexafluorid SF6

    // Computed pleural pressure
    double* pleuralPressure;
};

#endif

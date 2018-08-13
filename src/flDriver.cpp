#include "flDriver.h"
#include "lung.h"
#include "duct.h"
#include "acinus.h"
#include "gas.h"
#include "IOdict.h"
#include "dataStruct.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <eigen_3_3_4/Dense>
#include <eigen_3_3_4/Sparse>
#include <eigen_3_3_4/IterativeLinearSolvers>

using namespace std;
using namespace Eigen;



/// \file Driver function to run simulation
/// \param initData Contains all data of the lung input file
/// \param allbreathFlow Contains flow data of all breaths
/// \param controlProperties Control properties for simulation handle
/// \param systemProperties Properties of the fractal lung model and the trumpet acini
/// \param transportProperties  Transport properties (Diffusion coefficient, viscosity, density)
/// \return -
//**************************************************************/
void driver(init initData, breathFlow* allBreathFlow, breathResults* allBreathResults,
            controlProperties contProp, systemProperties sysProp, transportProperties transProp){

    // Initialize variables
    //**********************************************************/
    int k;

    double dt, Tend;
    double cumulativeTB;
    double Erel, Eabs, EQrel, EQabs;

    // Get variables from init data struct
    //**************************************************************/
    Tend = initData.Tend;
    dt = initData.dt;

    // Create Lung Object
    //**********************************************************/
    lung LUNG(&initData, &contProp, &sysProp, &transProp);

    // Generate morphology of conducting airways
    cout << "\n\nGenerate Morphology:\n";
    cout << "*******************************\n";
    LUNG.genMorphCondAirways(LUNG.pTrachea);

    // Compute total airway volume
    // (also set min/max generation of end duct and calculate total airway length)
    LUNG.airwayVolume(LUNG.pTrachea);
    cout<<std::left<<setw(20) << "# ducts" << " = " << LUNG.nbrDucts << endl;
    cout<<std::left<<setw(20) << "# end ducts" << " = " << LUNG.nbrEndDucts << endl;
    cout<<std::left<<setw(20) << "min/max generation" << " = " << LUNG.endDuctMinGen << "/" << LUNG.endDuctMaxGen << endl;

    // Size acinus template to match FRC and TAWV+TAV
    LUNG.sizeAcinusTemplate();

    // Generate properly scaled acini
    LUNG.genTrumpetAcinusOnAirways(LUNG.pTrachea);

    // Compute total acinus volume and acinus total resistance
    LUNG.acinusVolumeNResistance(LUNG.pTrachea);
    cout<<std::left<<setw(20) <<"FRC residual = " << abs(LUNG.TAWV+LUNG.TAV - LUNG.FRC) << endl;
    cout<<std::left<<setw(20) <<"TAV / FRC    = " << abs(LUNG.TAV)/abs(LUNG.TAWV+LUNG.TAV) << endl;

    cout << "\n\nRead files ..." << endl;
    // Read breath periods and tidal volume corresponding to measured flow signal
    LUNG.readTBTV(allBreathFlow);
    LUNG.setTBTV(0);

    // Read files for boundary conditions
    if (LUNG.bc==0){
        // Read generic data for inlet flow
        LUNG.readInletFlow(allBreathFlow);
    }
    else {
        // Read generic data for inlet flow
        LUNG.readPleuralPressure(allBreathFlow);
    }

    // Read table with transmissability correction factors due to pulsatile flow
    LUNG.readTransFact();

    // Read table with modification values and modify Lung morphology
    LUNG.readModifications();
    LUNG.applyModifications(LUNG.pTrachea);
    LUNG.correctTotalVolumeAfterModification(LUNG.pTrachea);
    cout << "\n\nModification table:"<< endl;
    cout << "*******************************" << endl;
    cout << LUNG.modTable << "\n" << endl;
    cout << "Lobular volume ratio pre / post modifications: " << LUNG.TAV/LUNG.TAV_aft_mod << endl;

    // Calculate trumpet acinus stiffness
	  LUNG.calculateStiffnessParametersInTrumpetAcinus(LUNG.pTrachea);

    // Initialize matrices and vectors for liner system for pressure
    LUNG.initializeMat();
    LUNG.airwaySetTK1(LUNG.pTrachea);
    LUNG.writeCoeffMatrix(LUNG.pTrachea);

    // Pass coefficient matrix to solver
    //cout << "\nDecomposing coefficient matrix ..." << endl;
    //ColPivHouseholderQR<MatrixXd> solver(LUNG.coeffMat);

     // Setup iterative solver
    Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double> > solverBiCGSTAB;
    solverBiCGSTAB.preconditioner().setDroptol(.001);

    // Estimate velocity distribution for grid size & (total) number of ducts
    cout << "\nEstimating velocity distribution at peak inlet flow rate ..." << endl;
    LUNG.writeRHSVector(LUNG.pTrachea, 1);
    // LUNG.pressureVec = solver.solve(LUNG.rhsVec); // used together with matrix decompositoin e.g. ColPivHousholder
    LUNG.pressureVec = solverBiCGSTAB.compute((LUNG.coeffMat).sparseView()).solve(LUNG.rhsVec);
    LUNG.distrPressure(LUNG.pTrachea);
    LUNG.estimateFlow(LUNG.pTrachea);

    // Setup species
    cout << "\n\nSetup gas ..." << endl;
    LUNG.setupGaseousSpecies(LUNG.pTrachea);
    LUNG.computeTotalNbrOfGridPoints(LUNG.pTrachea);

    // estimate maximum time-step refinement
    LUNG.computeTimeStepRefinement(LUNG.pTrachea, 0);
    cout << "maximum time step refinement = " << LUNG.Nt << endl;

    // Initialize time series
    LUNG.initializeVTKArrays();

    // Start time loop
    //**********************************************************/
    cout << "\n\nStart time loop ..." << endl;
    cout << "*******************************" << endl;
    cumulativeTB = k = 0;

    for (int it = 0; it<int(round(Tend/dt)); it++){

        // Check for new breath
        if (it==int(round(cumulativeTB/dt))){
            // Set new temporary TB and TV
            LUNG.setTBTV(k);
            cumulativeTB += LUNG.TB;
            k++;
        }

        // compute time step refinement for upcoming pressure, flow and species species concentration update (in order to satisfy CFL condition)
        LUNG.Nt = 1.0; //re-initialize
        LUNG.computeTimeStepRefinement(LUNG.pTrachea, 1);

        // Sub time stepping to update acini and ducts
        for (int n = 0; n<LUNG.Nt; n++){

            // Update entries in coefficient matrix
            LUNG.airwaySetTK1(LUNG.pTrachea);
            LUNG.coeffMat.setZero();
            LUNG.writeCoeffMatrix(LUNG.pTrachea);

            // Write right hand side
            LUNG.interpInputData(it, n);
            LUNG.writeRHSVector(LUNG.pTrachea, 0);

            // Solve for pressure
            // LUNG.pressureVec = (LUNG.coeffMat).colPivHouseholderQr().solve(LUNG.rhsVec);
            LUNG.pressureVec = solverBiCGSTAB.compute((LUNG.coeffMat).sparseView()).solve(LUNG.rhsVec);

            // Distribute pressure in airway tree
            LUNG.distrPressure(LUNG.pTrachea);
            if (LUNG.bc==0){
                LUNG.ppl = LUNG.pressureVec(LUNG.nbrDucts);
            }

            // Compute flow in each airway duct
            LUNG.Qouttot *= 0.;
            LUNG.computeFlow(LUNG.pTrachea);

            // Update flow properties in breathing acinus
            LUNG.updateAcinus(LUNG.pTrachea);

            // Solve transport equation for one time-step (Crank-Nicolson)
            LUNG.updateConcentrationInDucts(LUNG.pTrachea);

            // Update old concentration
            LUNG.updateOldConc(LUNG.pTrachea);
        }

        // Write error of solution for pressure problem and error of flow balance as well as time step refinement to terminal
        if (it%10 == 0){
            // Time step refinement
            cout << "time = " << it*dt << ", Time step refinement = " << LUNG.Nt << endl;

            // Pressure problem
            Erel = (LUNG.coeffMat*LUNG.pressureVec - LUNG.rhsVec).norm() / LUNG.rhsVec.norm();
            Eabs = (LUNG.coeffMat*LUNG.pressureVec - LUNG.rhsVec).norm();
            cout << "   "<< "Erel = " << Erel << "; Eabs = " << Eabs << endl;

            // Flow balance
            EQrel = (LUNG.Qin - LUNG.Qouttot)/LUNG.Qin;
            EQabs = (LUNG.Qin - LUNG.Qouttot);
            cout << "   " << "EQrel = " << EQrel << "; EQabs = " << EQabs << endl;
        }

        // Write Primary Output
        //******************************************************/
        LUNG.writePrimaryLungData(it, allBreathFlow, allBreathResults, LUNG.pTrachea);

        // Write full output data
        //******************************************************/
        if (contProp.writeFull && it%contProp.fOutFull == 0){
            LUNG.writeFullLungData(LUNG.pTrachea);
        }
    }
}

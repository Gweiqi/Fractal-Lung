// Include libraries and object files
//**************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "flDriver.h"
#include "dataStruct.h"
#include "IOdict.h"
#include <streambuf>
#include <ctime>

using namespace std;

int main(int argc, char** argv){

    double exe_time;
    clock_t start_clock, end_clock;

    // Print out info
    //**************************************************************/
    start_clock = clock();
    cout << "*******************************\n";
    cout << "...start simulation run\n\n";
    cout << "*******************************\n\n\n";


    // Initialize variables for read in
    //**********************************************************/
    int nbrBreaths = 0;
    double fs = 0;
    ifstream inFile;
    ofstream outFile;


    // Read initial data
    // (to be provided by data/inputLung)
    //**********************************************************/
    init initData;
    cout << "Input parameters:\n";
    cout << "*******************************\n";
    IOdict initDataIO    = IOdict("data/inputLung");
    initDataIO.outputConsole = 1;
    initData.dL         = initDataIO.lookup("dLimit");
    initData.FRC        = initDataIO.lookup("FRC");
    initData.washout    = initDataIO.lookup("washout");
    initData.species    = initDataIO.lookup("species");

    // Read lung parameter from external file
    // (to be provided by data/systemProperties)
    //**********************************************************/
    systemProperties sysProp;
    cout << "\n\nLung parameters:\n";
    cout << "*******************************\n";
    IOdict sysPropIO            = IOdict("constant/systemProperties");
    sysPropIO.outputConsole     = 1;
    sysProp.maxGenLb            = sysPropIO.lookup("maxGenLb");
    sysProp.p0                  = sysPropIO.lookup("p0");
    sysProp.lLbm_fac             = sysPropIO.lookup("lLbm_fac");
    sysProp.bending             = sysPropIO.lookup("bending");
    sysProp.r                   = sysPropIO.lookup("r");
    sysProp.eta                 = sysPropIO.lookup("eta");
    sysProp.kappa               = sysPropIO.lookup("kappa");
    sysProp.n1                  = sysPropIO.lookup("n1");
    sysProp.n2                  = sysPropIO.lookup("n2");
    sysProp.z_star              = sysPropIO.lookup("z_star");
    sysProp.cin                 = sysPropIO.lookup("cInlet");
    sysProp.bc                  = sysPropIO.lookup("bc");
    sysProp.R_fac               = sysPropIO.lookup("R_fac");
    sysProp.Diff_fac            = sysPropIO.lookup("Diff_fac");
    sysProp.scalingLbL          = sysPropIO.lookup("scalingLbL");


    // Read fluid parameter from external file
    // (to be provided by data/transportProperties)
    //**********************************************************/
    transportProperties transProp;
    cout << "\n\nFluid parameters:\n";
    cout << "*******************************\n";
    IOdict transPropIO          = IOdict("constant/transportProperties");
    transPropIO.outputConsole   = 1;
    transProp.mu                = transPropIO.lookup("mu");
    transProp.rho               = transPropIO.lookup("rho");
    transProp.Diff              = transPropIO.lookupVect("D_O2");


    // Read in number of breaths and sampling frequency
    // (to be provided by data/nbfs)
    //**********************************************************/
    inFile.open("data/nbfs");
    if (inFile.is_open()){
        inFile >> nbrBreaths;
        inFile >> fs;
    }
    inFile.close();

    // Create flow struct array and define it's members
    breathFlow* allBreathFlow = new breathFlow[nbrBreaths];


    // Read in breath periode (TB), tidal volume (TV) and number of data points (N)
    //**********************************************************/
    inFile.open("data/TBTVN");
    if (inFile.is_open()){
        for (int i = 0; i<nbrBreaths; i++){
            inFile >> allBreathFlow[i].TB;
            inFile >> allBreathFlow[i].TV;
            inFile >> allBreathFlow[i].N;
            // Create array for N data points
            allBreathFlow[i].flowData = new double[allBreathFlow[i].N];
            allBreathFlow[i].pplData  = new double[allBreathFlow[i].N];
            // Define sample frequency
            allBreathFlow[i].fs = fs;
        }
        inFile.close();
    }


    // Read in flow or pressure data
    //**********************************************************/
    if (sysProp.bc == 0){
        inFile.open("data/inletFlow");
        if (inFile.is_open()){
            for (int i = 0; i<nbrBreaths; i++){
                for (int j = 0; j<allBreathFlow[i].N; j++){
                    inFile >> allBreathFlow[i].flowData[j];
                }
            }
            inFile.close();
        }
    }
    else{
        inFile.open("data/pleuralPressure");
        if (inFile.is_open()){
            for (int i = 0; i<nbrBreaths; i++){
                for (int j = 0; j<allBreathFlow[i].N; j++){
                    inFile >> allBreathFlow[i].pplData[j];
                }
            }
            inFile.close();
        }
    }



    // Time integration step is the inverse of the sample frequency of the flow
    //**********************************************************/
    double dt, Tend;
    dt = 1./allBreathFlow[0].fs;
    initData.dt = dt;

    // Compute total simulation time;
    Tend = 0.;
    for (int i = 0; i<nbrBreaths; i++){
        Tend += allBreathFlow[i].TB;
    }
    initData.Tend = Tend;
    initData.nbrBreaths = nbrBreaths;


    // Read/asign control parameter from external file
    // (to be provided by data/controlProperties)
    //**********************************************************/
    controlProperties contProp;
    cout << "\n\nControl properties:\n";
    cout << "*******************************\n";
    IOdict contPropIO           = IOdict("constant/controlDict");
    contPropIO.outputConsole    = 1;
    contProp.dt                 = dt; cout<<std::left<<setw(20)<<"dt"<<" = "<<dt<<endl;
    contProp.Tend               = Tend; cout<<std::left<<setw(20)<<"Tend"<<" = "<<Tend<<endl;
    contProp.nbrBreaths         = nbrBreaths; cout<<std::left<<setw(20)<<"#Breaths"<<" = "<<nbrBreaths<<endl;
    contProp.generalCN          = contPropIO.lookup("generalCN");
    contProp.CFL                = contPropIO.lookup("CFL");
    contProp.NxLob              = contPropIO.lookup("NxLob");
    contProp.NxDuctMin          = contPropIO.lookup("NxDuctMin");
    contProp.NxDuctMax          = contPropIO.lookup("NxDuctMax");
    contProp.fOutFull           = contPropIO.lookup("fOutFull");
    contProp.writeFull          = contPropIO.lookup("writeFull");



    // Prepare results array for species concentration at inlet and pleural pressure
    //**********************************************************/
    breathResults* allBreathResults = new breathResults[nbrBreaths];
    for (int i = 0; i<nbrBreaths; i++){
        allBreathResults[i].N = allBreathFlow[i].N;

        // Create arrays for N data points
        allBreathResults[i].time  = new double[allBreathResults[i].N];
        allBreathResults[i].speciesIData  = new double[allBreathResults[i].N];
        allBreathResults[i].speciesIIData = new double[allBreathResults[i].N];
        allBreathResults[i].pleuralPressure = new double[allBreathResults[i].N];
    }

    // Call driver for simulation start
    //**********************************************************/
    driver(initData, allBreathFlow, allBreathResults, contProp, sysProp, transProp);


    // Write simulated data to file (primary output)
    //**********************************************************/
    remove("data/primary_results");
    outFile.open("data/primary_results");
    outFile.precision(10);
    outFile.setf(std::ios::fixed);

    for (int i = 0; i<nbrBreaths; i++){
        for (int j = 0; j<allBreathResults[i].N; j++){
            if (allBreathResults[i].time[j] < 17){
                //cout << "i = " << i << "; j = " << j << "; t = " << allBreathResults[i].time[j] << "; c = " << allBreathResults[i].speciesIData[j] << "; ppl = " << allBreathResults[i].pleuralPressure[j] << endl;
            }
            outFile<<allBreathResults[i].time[j]<<"\t"<< allBreathResults[i].speciesIData[j]<<"\t"<< allBreathResults[i].speciesIIData[j]<<"\t"<< allBreathResults[i].pleuralPressure[j]<<"\n";
        }
    }
    outFile.close();

    //  Clean up variables
    //**********************************************************/
    delete[] allBreathFlow;
    delete[] allBreathResults;

    end_clock = clock();
    exe_time = float(double(end_clock - start_clock)/CLOCKS_PER_SEC);
    cout << "\nDone" << endl;
    cout << "Execution time: " << exe_time << " s (" << exe_time/60. << " min, " << exe_time/3600. << " h)" << endl;

    return 0;
}

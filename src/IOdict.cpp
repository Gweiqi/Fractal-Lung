#include "IOdict.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <iomanip>


// Constructor
//**************************************************************/
IOdict::IOdict(){
    fileName = "";
    outputConsole = 0;
}

IOdict::IOdict(string var){
    fileName = var;
    outputConsole = 0;
}


// Destructor
//**************************************************************/
IOdict::~IOdict(){
}


// Member Functions
//**************************************************************/

// Print current file name
void IOdict::printFileName(){

}


// Clear file
void IOdict::clearFile(){
    remove(fileName.c_str());
}


// Search for specific variable, return value
double IOdict::lookup(string varName){
    string line, lineWord;
    ifstream file(fileName.c_str());

    if(file.is_open()){

        while(file.good()){
            getline(file,line);
            std::stringstream lineStream(line);

            // Go throughout whole line
            while (lineStream >> lineWord){
                // Find Begin of comment and break
                if(lineWord.substr(0,2)=="//"){
                    break;
                }

                // Check for searched variable name, return it
                if(varName==lineWord){
                    lineStream >> lineWord;
                    if(outputConsole){
                        cout<<std::left<<setw(20)<<varName<<" = "<<lineWord<<endl;
                    }
                    file.close();
                    return atof(lineWord.c_str());
                }
            }
        }
        // In case variable was not found:
        cout<<"Variable '"<< varName<<"' not found"<<endl;
        file.close();
        return 0;
    }

    // In case variable was not found:
    cout<<"File '"<< fileName<<"' not found"<<endl;
    file.close();
    return 0;
}


// Search for specific vector, return vector
vector<double> IOdict::lookupVect(string varName){
    int i = 0;
    string line, lineWord;
    string test;
    ifstream file(fileName.c_str());
    vector<double> vect;

    if(file.is_open()){

        while(file.good()){
            getline(file,line);
            std::stringstream lineStream(line);

            // Go throughout whole line
            while (lineStream >> lineWord){
                // Find Begin of comment and break
                if(lineWord.substr(0,2)=="//"){
                    break;
                }

                // Check for searched variable name, return it
                if(varName==lineWord){
                    // Go through line by word
                    while(lineStream >> lineWord){
                        // Break when comment reached
                        if(lineWord.substr(0,2)=="//"){
                            break;
                        }
                        // Otherwise store data to vector
                        else{
                        vect.push_back(atof(lineWord.c_str()));
                        i++;
                        }
                    }

                    // Output to console
                    if(outputConsole){
                         cout<<std::left<<setw(20)<<varName<<" = [";
                        for (int i=0; i<vect.size();i++){
                            cout << vect[i]<<" ";
                        }
                        cout<<"]"<<endl;
                    }
                    file.close();
                    return vect;
                }
            }
        }

        // In case variable was not found:
        cout<<"Variable '"<< varName<<"' not found"<<endl;
        file.close();
        return vect;
    }

    // In case variable was not found:
    cout<<"File '"<< fileName<<"' not found"<<endl;
    file.close();
    return vect;
}



// Write vector to file
void IOdict::write(VectorXd data){
    std::ofstream file(fileName.c_str(),std::ofstream::app);
    if(file.is_open()){
        file<<data.transpose()<<endl;
    }
    file.close();
}


// Write array to file
//**************************************************************/
void IOdict::write(double *data, int length){
    std::ofstream file(fileName.c_str(),std::ofstream::app);
    file.precision(8);
    file.setf(std::ios::fixed);

    if(file.is_open()){
        for(int i=0;i<length;i++){
            file<<data[i]<<" ";
        }
        file<<endl;
    }
    file.close();
}


// Write data to file (1 specie)
//**************************************************************/
void IOdict::write(double data){
    std::ofstream file(fileName.c_str(),std::ofstream::app);
    file.precision(8);
    file.setf(std::ios::fixed);

    if(file.is_open()){
        file<<data<<endl;
    }
    file.close();
}

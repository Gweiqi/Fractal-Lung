#ifndef IODICT
#define IODICT

#include <iostream>
#include <vector>
#include <eigen_3_3_4/Dense>

using namespace std;
using namespace Eigen;

/// IOdict
class IOdict {
public:
    // Public variables
    //**********************************************************/
    /// Bool to print on console
    bool outputConsole;

    /// Filename
    string fileName;

    // Public member functions
    //**********************************************************/
    /// Constructor
    IOdict(string fileName);
    IOdict();

    /// Destructor
    ~IOdict();

    /// Print current file name
    void printFileName();

    /// Clear file
    void clearFile();

    /// Search for specific variable, return value
    double lookup(string varName);

    /// Search for specific vector, return vector
    vector<double> lookupVect(string varName);

    /// Write vector into file
    void write(VectorXd data);

    /// Write array into file
    void write(double *data, int length);

    /// Write double into file
    void write(double data);
};


#endif

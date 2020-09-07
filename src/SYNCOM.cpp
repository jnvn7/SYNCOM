// SYNCOM - A Nonlinear Synthetic Rope Numerical Computation Software
//
// Copyright 2020 Jessica Nguyen <nvnguyen@umass.edu>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
////////////////////////////////////////////////////////////////////////////////

#include "error.h"
#include "readIn.h"
#include "setting.h"
#include "strainSolver.h"
#include "stressSolver.h"
#include "printOut.h"
#include <chrono>
#include <string>
#include <iostream>

using namespace std;
using namespace std::chrono;
using namespace rope;

////////////////////////////////////////////////////////////////////////////////
/// Check the generated errors. Close the application window if an error is catched.
////////////////////////////////////////////////////////////////////////////////
void catch_error(ErrorCode errCodes, ErrorOut errOut)
{
    if (errCodes == ErrorCode::SUCCESS || errCodes == ErrorCode::SIMULATION_COMPLETED) {
        cout << "    ... " + errOut.message(errCodes) << endl << endl;
    } 
    else {
        cout << "  ...Error: " + errOut.message(errCodes) << endl << endl;

        char usr_input[2];
        printf("\n\n Exit Application? (y) \n\n");
        char s[2];
#ifndef __unix__
        scanf_s("%s", usr_input, (unsigned)_countof(s));
#else
        scanf("%s", usr_input);
#endif
        if (strcmp(usr_input, "y") == 0) exit(1);
    }   
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the simulation.
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    /// Starting time;
    auto start = high_resolution_clock::now();

    /// Print copyright to screen;
    print_copyright();

    /// Initialize instances for read, write and errors;
    cout << "Initialization..." << endl << endl;
    ErrorCode errCodes;
    ErrorOut errOut;
    ReadIn readInput;

    // Set work folder and input files.
    MatProps mat_props;     // stores material properties.
    Setting setting(argv[0], &mat_props);

    /// Initialization of setting from Setting.xml.
    cout << "Reading Setting.xml ..." << endl;
    errCodes = readInput.readIn_data(setting);
    catch_error(errCodes, errOut);
    cout << "Validating input data..." << endl;
    errCodes = setting.validate();
    catch_error(errCodes, errOut);
  
    /// Starts simulation;
    cout << "Simulation starts..." << endl;

    if (setting.module == 0) {
        // Initialization of output instance;
        strainSolver strainOutput(setting);

        // Perform the simulation
        errCodes = strainOutput.syncom_solver(setting);

        // Write results to file;
        print_mod1(strainOutput, setting);
    }
    else {
        // Initialization of output instance;
        stressSolver stressOutput(setting);
        
        // Perform the simulation
        errCodes = stressOutput.syncom_solver(setting);

        // Write results to file;
        print_mod2(stressOutput, setting);
    }
    
    catch_error(errCodes, errOut);
    
    /// Stopping time and duration;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Computation completed in: " << duration.count() << "(milliseconds) " 
        << "for " << setting.dataIn.size() << " time steps!" << endl;
  
    end_simu();
} 

// End of main program

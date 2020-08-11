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
/// Print description on the screen.
////////////////////////////////////////////////////////////////////////////////
void print_copyright(void)
{
    cout << endl << endl;
    cout << "------------------------------------------------------------------"
        "--------------" << endl << endl;
    cout << "           A Nonlinear Synthetic Constitutive Modeling Tool\n";
    fprintf(stdout, " %s Version 1.0     \n", "                         SYNCOM");
    cout << "                  Copyright (c) 2020 Jessica Nguyen\n\n";
    cout << "------------------------------------------------------------------"
        "--------------" << endl << endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the generated erros. Quit if an error is found.
////////////////////////////////////////////////////////////////////////////////
void check_status(ErrorCode errCodes, ErrorOut errOut)
{
    if (errCodes == ErrorCode::SUCCESS || errCodes == ErrorCode::SIMULATION_COMPLETED) {
        cout << "    ... " + errOut.message(errCodes) << endl << endl;
    } 
    else {
        cout << "  ... Unseccessful: " + errOut.message(errCodes) << endl << endl;

        char usr_input[2];
        printf("\n\n Exit Application? (y/n) \n\n");
        char s[2];
        scanf_s("%s", usr_input, (unsigned)_countof(s));

        if (strcmp(usr_input, "y") == 0) exit(1);
    }   
}

////////////////////////////////////////////////////////////////////////////////
/// End of Simulation.
////////////////////////////////////////////////////////////////////////////////
void end_simu()
{
    std::ios_base::sync_with_stdio(false);
    std::cout << "\nPress enter to close the application...\n";
    std::cin.ignore();
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
    cout << "  Reading Setting.xml ..." << endl;
    errCodes = readInput.readIn_data(setting);
    check_status(errCodes, errOut);
    cout << "  Validating input data..." << endl;
    errCodes = setting.validate();
    check_status(errCodes, errOut);
  
    /// Starts simulation;
    cout << "  Simulation starts..." << endl;

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

    check_status(errCodes, errOut);

    /// Stopping time and duration;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Computation completed in: " << duration.count() << "(milliseconds) " 
        << "for " << setting.dataIn.size() << " time steps!" << endl;
  
    end_simu();
} 

// End of main program

// SYNCOM - A Nonlinear Synthetic Rope Numerical Computation Software
//
// Copyright (c) 2020 Jessica Nguyen <nvnguyen@umass.edu>
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

#include "SynCOM.h"

//using namespace std;

namespace rope {

    ////////////////////////////////////////////////////////////////////////////////
    /// Read User's Inputs from MaterDef.xml and initilize the application.
    ////////////////////////////////////////////////////////////////////////////////
    int SC_initialize(int numNodes, string input_file, string line_type, stressSolver& stressSolver) {

        /// Sets up working folder;
        ReadIn readInput;
        ErrorCode errCodes;
        ErrorOut errorOut;
        std::string filename(input_file);
        stressSolver.sc_folder = filename.substr(0, filename.find_last_of("/\\") + 1);
        stressSolver.sc_inputfile = filename;
        
        /// Print info and read inputs;
        errCodes = readInput.readIn_data(line_type, stressSolver);

        // Print copy_right;
        print_copyright(stressSolver);
        print_message(stressSolver, "  Reading MaterDef.xml ...\n");

        // Read inputs;
        if (errCodes != ErrorCode::SUCCESS)
        {
            print_log(stressSolver, errCodes, errorOut);
            return 1;
        }
        else {
            print_message(stressSolver, "   .....Successful\n\n");
        }

        print_message(stressSolver, "  Validating input data...\n");
        errCodes = stressSolver.validate();
        if (errCodes != ErrorCode::SUCCESS)
        {
            print_log(stressSolver, errCodes, errorOut);
            return 1;
        }
        else {
            print_message(stressSolver, "   .....Successful\n\n");
        }
        
        // Initializes values of object from input data;
        stressSolver.stressSolver_init(numNodes);
        
        return 0;

    } // End of initialize func.

    ////////////////////////////////////////////////////////////////////////////////
    /// End simulation.
    ////////////////////////////////////////////////////////////////////////////////
    int close_app(stressSolver& stressSolver) {
        print_message(stressSolver, "   .....Computation completed!\n");
        return 0;

    } // End of close_app

    ////////////////////////////////////////////////////////////////////////////////
    /// Print log file for error check.
    ////////////////////////////////////////////////////////////////////////////////
    int print_log(stressSolver& stressSolver, ErrorCode errCodes, ErrorOut errOut)
    {
#ifndef __unix__
        char buffer[32];
#if _WIN64
        _time64(&stressSolver.aclock);   // Get time in seconds.
        _localtime64_s(&stressSolver.newtime, &stressSolver.aclock);   // Convert time to struct tm form.
#else
        _time32(&stressSolver.aclock);   // Get time in seconds.
        _localtime32_s(&stressSolver.newtime, &stressSolver.aclock);   // Convert time to struct tm form.
#endif
        asctime_s(buffer, 32, &stressSolver.newtime);
#endif 

        std::ofstream log_file;
        log_file.open(stressSolver.log_filename, std::ofstream::app);

#ifndef __unix__
        log_file << "  " << buffer << errOut.message(errCodes) << endl << endl;
#else
        log_file << "  " << errOut.message(errCodes) << endl << endl;
#endif
        log_file.close();
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////7
    /// Write log file for error check.
    ////////////////////////////////////////////////////////////////////////////////
    int print_message(stressSolver& stressSolver, string message)
    {
#ifndef __unix__
        char buffer[32];
#if _WIN64
        _time64(&stressSolver.aclock);   // Get time in seconds.
        _localtime64_s(&stressSolver.newtime, &stressSolver.aclock);   // Convert time to struct tm form.
#else
        _time32(&stressSolver.aclock);   // Get time in seconds.
        _localtime32_s(&stressSolver.newtime, &stressSolver.aclock);   // Convert time to struct tm form.
#endif
        asctime_s(buffer, 32, &stressSolver.newtime);
#endif 

        std::ofstream log_file;
        log_file.open(stressSolver.log_filename, std::ofstream::app);

#ifndef __unix__
        log_file << buffer << message << endl << endl;
#else
        log_file << message << endl << endl;
#endif
        log_file.close();
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Print copyright.
    ////////////////////////////////////////////////////////////////////////////////
    void print_copyright(stressSolver& stressSolver)
    {
        string message = "---------------------------------------------------------"
            "------------\n"
            "           A Nonlinear Synthetic Constitutive Modeling Tool \n"
            "                         SYNCOM Version 1.0                 \n"
            "                  Copyright (c) 2020 Jessica Nguyen         \n"
            "---------------------------------------------------------------------\n";

        ofstream log_file;
        log_file.open(stressSolver.log_filename, std::ofstream::trunc);
        log_file << message << endl << endl;
        log_file.close();
    }

} // End of namespace rope

// End of main program
  
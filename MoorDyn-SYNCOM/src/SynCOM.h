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


#ifndef SynCOM_h
#define SynCOM_h

#include "SC_error.h"
#include "SC_readIn_api.h"
#include "SC_stressSolver_api.h"
#include <string>
#include <iostream>
#include <time.h>
#include <fstream>
#include <stdio.h>

namespace rope {
    extern 

    // Initialization;
    int SC_initialize(int numNodes, string input_file, string line_type, stressSolver& stressSolver);
   
    // Simulation completed;
    int close_app(stressSolver& stressSolver);

    // Print log file for error check;
    int print_log(stressSolver& stressSolver, ErrorCode errCodes, ErrorOut errOut);

    // Write log file for error check;
    int print_message(stressSolver& stressSolver, string message);

    // Print copyright;
    void print_copyright(stressSolver& stressSolver);

/*
    // Clear all global variables and close the program.
    int finish(void); */


}

#endif // SYNCOM_h


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


#ifndef printOut_h
#define printOut_h

#include "strainSolver_api.h"
#include "stressSolver_api.h"
#include "error.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>

namespace rope {

    /// \brief Unitility functions.
    int print_message(string& filename, string message);
    void print_copyright(string& filename);

    /// NEED FIXING - Write log file.
    int print_log(string& file_name, ErrorCode errCodes, ErrorOut errOut);

} // End of namespace rope.

#endif // writeOut_h
#pragma once

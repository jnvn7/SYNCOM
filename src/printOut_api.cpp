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

#include "printOut_api.h"

namespace rope {
    ////////////////////////////////////////////////////////////////////////////////
    /// Print log file for error check.
    ////////////////////////////////////////////////////////////////////////////////
    int print_log(string& filename, ErrorCode errCodes, ErrorOut errOut)
    {
        std::ofstream log_file;
        log_file.open(filename, std::ofstream::app);
        log_file << "  " << errOut.message(errCodes) << endl << endl;
        log_file.close();
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////7
    /// Write log file for error check.
    ////////////////////////////////////////////////////////////////////////////////
    int print_message(string& filename, string message)
    {
        std::ofstream log_file;
        log_file.open(filename, std::ofstream::app);
        log_file << message << endl << endl;
        log_file.close();
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Print copyright.
    ////////////////////////////////////////////////////////////////////////////////
    void print_copyright(string& filename)
    {
        string message = "---------------------------------------------------------"
            "--------------\n"
            "           A Nonlinear Synthetic Constitutive Modeling Tool\n"
            "                         SYNCOM Version 1.0                 \n"
            "                  Copyright (c) 2020 Jessica Nguyen"
            "------------------------------------------------------------------"
            "--------------";

        ofstream log_file;
        log_file.open(filename, std::ofstream::app);
        log_file << message << endl << endl;
        log_file.close();
    }

} // End of namespace rope.

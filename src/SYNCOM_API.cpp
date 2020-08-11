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
#include "readIn_api.h"
#include "setting.h"
#include "strainSolver_api.h"
#include "stressSolver_api.h"
#include "printOut_api.h"
#include "SYNCOM_API.h"
#include <chrono>
#include <string>
#include <iostream>

using namespace std;
using namespace std::chrono;
using namespace rope;

////////////////////////////////////////////////////////////////////////////////
/// Construct objects to handle data.
////////////////////////////////////////////////////////////////////////////////
ErrorCode errCodes = ErrorCode::SUCCESS;
ErrorOut errorOut;
MatProps mat_props;
Setting setting(&mat_props);
strainSolver strainOutput(setting);
stressSolver stressOutput(setting);

////////////////////////////////////////////////////////////////////////////////
/// Read User's Inputs from Setting.xml and initilize the application.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR initialize(int module, char input_file[]) {

    /// Sets up working folder;
    ReadIn readInput;
    std::string filename(input_file);
    setting.setting_folder = filename.substr(0, filename.find_last_of("/\\") + 1);
    setting.setting_file = filename;
    setting.module = module;

    /// Print info and read inputs;
    errCodes = readInput.readIn_data(setting);

    // Print copy_right;
    print_copyright(setting.log_filename);
    print_message(setting.log_filename, "  Reading Setting.xml ...\n");

    // Read inputs;
    if (errCodes != ErrorCode::SUCCESS)
    {
        print_log(setting.log_filename, errCodes, errorOut);
        return 1;
    }
    else {
        print_message(setting.log_filename, "   .....Successful\n\n");
    }

    print_message(setting.log_filename, "  Validating input data...\n");
    errCodes = setting.validate();
    if (errCodes != ErrorCode::SUCCESS)
    {
        print_log(setting.log_filename, errCodes, errorOut);
        return 1;
    }
    else {
        print_message(setting.log_filename, "   .....Successful\n\n");
    }

    // Initialization of output instance;
    if (module == 0)
        strainOutput = strainSolver(setting);
    else {
        stressOutput = stressSolver(setting);
    }

    /// Starts simulation;
    print_message(setting.log_filename, "  Simulation can start...\n");
    return 0;
 
} // End of initialize func.

////////////////////////////////////////////////////////////////////////////////
/// Perform the simulation.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR SYNCOM(double dataIn)
{
    if (setting.module == 0) 
        // Perform the simulation
        errCodes = strainOutput.syncom_solver(setting, dataIn);
    else 
        // Perform the simulation
        errCodes = stressOutput.syncom_solver(setting, dataIn);

    /// Check simulation end status.
    if (errCodes != ErrorCode::SUCCESS)
    {
        print_log(setting.log_filename, errCodes, errorOut);
        return 1;
    }

    return 0;

} // End of SYNCOM func.

////////////////////////////////////////////////////////////////////////////////
/// Extracts Stress || Strain (deformation) results.
////////////////////////////////////////////////////////////////////////////////
/// Total Strain;
double DECLDIR extract_eps(void)
{
    if (setting.module == 0)
        return strainOutput.get_eps();
    else
        return stressOutput.get_eps();
} 

///Visco-elastic strain;
double DECLDIR extract_eps_ve(void)
{
    if (setting.module == 0)
        return strainOutput.get_eps_ve();
    else
        return stressOutput.get_eps_ve();
}

/// Visco-plastic strain;
/// Total Strain;
double DECLDIR extract_eps_vp(void)
{
    if (setting.module == 0)
        return strainOutput.get_eps_vp();
    else
        return stressOutput.get_eps_vp();
}

/// Stress;
double DECLDIR extract_sigma(void)
{
    if (setting.module == 0)
        return strainOutput.get_sigma();
    else
        return stressOutput.get_sigma();
}

// End of main program

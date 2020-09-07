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

#define EXPORT_SYNCOM_API

#include "error.h"
#include "readIn_api.h"
#include "setting.h"
#include "strainSolver_api.h"
#include "stressSolver_api.h"
#include "printOut_api.h"
#include "SynCOM_API.h"
#include <chrono>
#include <string>
#include <iostream>

using namespace std;
using namespace std::chrono;
//using namespace rope;

////////////////////////////////////////////////////////////////////////////////
/// Construct objects to handle data.
////////////////////////////////////////////////////////////////////////////////
rope::ErrorCode errCodes = rope::ErrorCode::SUCCESS;
rope::ErrorOut errorOut;
rope::MatProps mat_props;
rope::Setting setting(&mat_props);
rope::strainSolver strainOutput(setting);
rope::stressSolver stressOutput(setting);

////////////////////////////////////////////////////////////////////////////////
/// Read User's Inputs from Setting.xml and initilize the application.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR initializeSC(int module, char input_file[]) {

    /// Sets up working folder;
    rope::ReadIn readInput;
    std::string filename(input_file);
    setting.setting_folder = filename.substr(0, filename.find_last_of("/\\") + 1);
    setting.setting_file = filename;
    setting.module = module;

    /// Print info and read inputs;
    errCodes = readInput.readIn_data(setting);
    
    // Print copy_right;
    print_copyright(setting);
    print_message(setting, "  Reading Setting.xml ...\n");

    // Read inputs;
    if (errCodes != rope::ErrorCode::SUCCESS)
    {
        print_log(setting, errCodes, errorOut);
        return 1;
    }
    else {
        print_message(setting, "   .....Successful\n\n");
    }

    print_message(setting, "  Validating input data...\n");
    errCodes = setting.validate();
    if (errCodes != rope::ErrorCode::SUCCESS)
    {
        print_log(setting, errCodes, errorOut);
        return 1;
    }
    else {
        print_message(setting, "   .....Successful\n\n");
    }

    // Initialization of output instance;
    if (module == 0)
        strainOutput = rope::strainSolver(setting);
    else {
        stressOutput = rope::stressSolver(setting);
    }

    /// Starts simulation;
    print_message(setting, "  Computation starts...\n");
    return 0;
 
} // End of initialize func.

////////////////////////////////////////////////////////////////////////////////
/// Perform the simulation.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR SynCOM(double dataIn, double dt)
{
    if (setting.module == 0) 
        // Perform the simulation
        errCodes = strainOutput.syncom_solver(setting, dataIn, dt);
    else 
        // Perform the simulation
        errCodes = stressOutput.syncom_solver(setting, dataIn, dt);

    /// Check simulation end status.
    if (errCodes != rope::ErrorCode::SUCCESS)
    {
        print_log(setting, errCodes, errorOut);
        return 1;
    }

    return 0;

} // End of SYNCOM func.

////////////////////////////////////////////////////////////////////////////////
/// End simulation.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR close_app(void) {
    print_message(setting, "   .....Computation completed!\n");
    return 0;
}

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

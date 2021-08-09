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


#include "SC_error.h"

namespace rope {

    std::string ErrorOut::message(ErrorCode err)
    {
        switch (err)
        {
        case ErrorCode::SUCCESS:
            return "Successful.";

        case ErrorCode::MATERDEF_FILE_NO_MATERDEF_NODE:
            return "Root node of MaterDef.xml must be MaterDef.";

        case ErrorCode::MATERDEF_FILE_ERROR_PARSE:
            return "Found xml parse error in MaterDef.xml.";

        case ErrorCode::MATERDEF_FILE_NONEXISTENT:
            return ("Cant find MaterDef.xml file.");

        case ErrorCode::MATERDEF_FILE_NO_MODULE_SELECTION:
            return ("No module specified in MaterDef.xml file.");

        case ErrorCode::MATERDEF_FILE_BAD_MODULE_SELECTION:
            return ("Bad module specified in MaterDef.xml file.");

        case ErrorCode::BAD_MATERIAL_PROPERTIES_INPUT:
            return ("Check material property inputs.");

        case ErrorCode::BAD_NUMERICAL_MATERDEF_INPUT:
            return ("Check numerical MaterDef input.");

        case ErrorCode::BAD_DT_INPUT:
            return ("Negative dt detected.");

        case ErrorCode::NEGATIVE_STRAIN_DETECTED:
            return ("Negative strain detected.");

        case ErrorCode::MATERDEF_FILE_NO_INPUT_DATA_FILE:
            return ("Cant find input data (stress or strain) file.");

        case ErrorCode::INPUT_DATA_FILE_NONEXISTENT:
            return ("Input data (stress or strain) file is not provided.");

        case ErrorCode::MATERDEF_FILE_MATERIAL_NAME_PROPERTIES:
            return ("Cant find Input (stress or strain) data file. Also check for material properties's name tag.");

        case ErrorCode::MATERDEF_FILE_INCOMPLETE_MATERIAL_PROPERTIES:
            return ("Insufficient material property input.");

        case ErrorCode::MATERDEF_FILE_NAN_MATERIAL_PROPERTIES:
            return ("Bad material property input.");

        case ErrorCode::MATERDEF_FILE_BAD_DN_VALUES:
            return ("Bad Dn's values input.");

        case ErrorCode::FAIL_TO_OPEN_INPUT_FILE:
            return ("Fail to open input file. Check input for data file.");

        case ErrorCode::WRONG_INPUT_FILE_FORMAT:
            return ("Only maximum of 1 columns are expected in the data file.");

        case ErrorCode::NON_LOGICAL_COEFFICIENT_INPUT:
            return ("Check input coefficients' step function setup.");

        case ErrorCode::SYNCOM_INITIALIZATION_FAILED:
            return ("SynCOM failed to find converged solution.");

        case ErrorCode::NAN_INPUT_DATA:
            return ("NaN values found in input data.");

        case ErrorCode::NEGATIVE_STRESS_INPUT:
            return ("Negative stress input - No compression allowed.");

        case ErrorCode::NEGATIVE_STRAIN_INPUT:
            return ("Negative strain input - No compression allowed.");

        case ErrorCode::NAN_OUTPUT:
            return ("Solver encounters NaN output.");

        case ErrorCode::NAN_OUTPUT_VISCO_ELASTIC_MODEL:
            return ("Visco-elastic Solver encounters NaN output.");

        case ErrorCode::NAN_OUTPUT_VISCO_ELASTIC_PLASTIC_MODEL:
            return ("Visco-elastic and Visco-plastic Solver encounters NaN output.");

        case ErrorCode::NO_CONVERGED_SOLUTION:
            return ("Can't find converged solution.");

        case ErrorCode::NO_CONVERGED_SOLUTION_VISCO_PLASTIC_SOLVER:
            return ("Can't find converged solution - Viscoplastic Solver.");

        case ErrorCode::SIMULATION_COMPLETED:
            return ("Simulation complete.");
        }
    }

} // End of namespace rope.

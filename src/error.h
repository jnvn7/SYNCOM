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


#ifndef error_h
#define error_h

#include <string>

namespace rope {

    enum class ErrorCode
    {
        /// Setting file's format and path.
        SUCCESS = 0,
        SETTING_FILE_NONEXISTENT,
        SETTING_FILE_ERROR_PARSE,
        SETTING_FILE_NO_SETTING_NODE,
        SETTING_FILE_NO_MODULE_SELECTION,
        SETTING_FILE_BAD_MODULE_SELECTION,

        ///  Check Material Properties' Inputs in Setting.h.
        BAD_MATERIAL_PROPERTIES_INPUT,
        BAD_NUMERICAL_SETTING_INPUT,
        BAD_DT_INPUT,
        FAIL_TO_OPEN_INPUT_FILE,
        WRONG_INPUT_FILE_FORMAT,
        NAN_INPUT_DATA,

        /// Check Material Properties' Inputs in ReadIn.h.
        SETTING_FILE_NO_INPUT_DATA_FILE,
        INPUT_DATA_FILE_NONEXISTENT,
        SETTING_FILE_NO_MATERIAL_PROPERTIES,
        SETTING_FILE_INCOMPLETE_MATERIAL_PROPERTIES,
        SETTING_FILE_NAN_MATERIAL_PROPERTIES,
        SETTING_FILE_BAD_DN_VALUES,

        /// Computation;
        NAN_OUTPUT,
        NAN_OUTPUT_VISCO_ELASTIC_MODEL,
        NAN_OUTPUT_VISCO_ELASTIC_PLASTIC_MODEL,
        NO_CONVERGED_SOLUTION,
        SIMULATION_COMPLETED
    };

    class ErrorOut
    {
    public:
        std::string message(ErrorCode err);
    };

} // End of namespace rope.

#endif //error_h

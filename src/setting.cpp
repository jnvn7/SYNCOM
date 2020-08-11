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


#include "setting.h"

namespace rope {

    Setting::Setting(const  std::string path, MatProps* mat_props) :
        app_path(path),
        setting_file("Setting.xml")
    {
        size_t folder_index = app_path.find_last_of("/\\");
        setting_folder = app_path.substr(0, folder_index + 1);

        material_props = mat_props;
    }


    ErrorCode Setting::validate(void)
    {
        if (material_props->sigma_yield0 < 0
            || material_props->MBL <= 0
            || material_props->Do <= 0) {
            return ErrorCode::BAD_MATERIAL_PROPERTIES_INPUT;
        }
        else {
            return ErrorCode::SUCCESS;
        }

        if (dt <= 0
            || limit <= 0
            || tol < 0) {
            return ErrorCode::BAD_NUMERICAL_SETTING_INPUT;
        }
        else {
            return ErrorCode::SUCCESS;
        }
    }

} // End of namespace rope.

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


#ifndef setting_h
#define setting_h

#include "error.h"
#include <vector>
#include <sstream>

namespace rope {

    struct MatProps
    {
        /// Material properties 
        // (usually obtained from experimental data);
        double   sigma_yield0;
        double   MBL;
        double   Do;
        double   sumDn = 0;
        std::vector<double> Dn;
        std::vector<double> lamdaN;
        std::vector<double> a0Coefs;
        std::vector<double> g0Coefs;
        std::vector<double> g1Coefs;
        std::vector<double> g2Coefs;
        std::vector<double> EpCoefs;
        std::vector<double> npCoefs;
        std::vector<double> H_vpCoefs;

        /// Stores coefficient step limit;
        std::vector<std::vector<double>> a0stress_lim;
        std::vector<std::vector<double>> g0stress_lim;
        std::vector<std::vector<double>> g1stress_lim;
        std::vector<std::vector<double>> g2stress_lim;
        std::vector<std::vector<double>> npstress_lim;
        std::vector<std::vector<double>> Epstress_lim;
        std::vector<std::vector<double>> Hstress_lim;

        std::vector<int> step_num;
    };

    /// \brief Setting sets up solver.
    ///
    /// Setting setups file arrangement and analysis options.
    class Setting
    {

    public:

        Setting(const  std::string path, MatProps* mat_props);
        Setting(MatProps* mat_props) {
            material_props = mat_props;
            material_props->step_num = std::vector<int>(7, 0);
        };

        ErrorCode validate(void);

        /// The full path of the compiled application SYNCOM.exe or SYNCOM.app.
        const  std::string app_path;

        /// Important: The Setting.xml should be in the same folder as the Applicatiom
        /// SYNCOM.
        //const  std::string setting_file;
        std::string setting_file;

        /// The folder is where the Setting.xml located and the other input files
        /// will be find in the subfolders and the output will be written into
        /// files located in the subfolders as well.
        std::string setting_folder;

        /// Extracted from full path for the input file.
        std::string work_folder;

        MatProps* material_props;

        /// Path to input file for stress/strain;
        std::string input_data_path;

        /// Path to output file for results;
        std::string output_filename;

        /// Path to log file and time parameters;
        std::string log_filename;

#ifndef __unix__
        struct tm newtime;
        __time32_t aclock;
#endif

        /// Module selection, time and Stress/Strain input data;
        int limit, module;
        double tol;
        std::vector<double> dt;
        std::vector<std::vector<double>> dataIn;
    };

} // End of namespace rope.

#endif // setting_h
#pragma once

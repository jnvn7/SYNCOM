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

#include "printOut.h"

namespace rope {
    ////////////////////////////////////////////////////////////////////////////////
    // Print results to file for Mod 1
    ////////////////////////////////////////////////////////////////////////////////
    void print_mod1(strainSolver& strainSolver, Setting& setting) {

        FILE* output_file;
        string name_ext = "_mod1.csv";
        fopen_s(&output_file, (setting.output_filename + name_ext).c_str(), "w");
        fprintf(output_file, "Time(s)       Stress(in)        Total_Strain(out)    "
                "   Visco-elastic_Strain    Visco-plastic_Strain \n");

        // Write cable state.
        for (size_t i = 0; i < strainSolver.eps.size(); i++)
        {
            fprintf(output_file, "% .5E % .5E % .5E % .5E % .5E \n",
                strainSolver.simTime[i], setting.dataIn[i][0], strainSolver.eps[i],
                strainSolver.eps_ve[i], strainSolver.eps_vp[i]);
        }
        fclose(output_file);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Print results to file for Mod 2
    ////////////////////////////////////////////////////////////////////////////////
    void print_mod2(stressSolver& stressSolver, Setting& setting) {

        FILE* output_file;
        string name_ext = "_mod2.csv";
        fopen_s(&output_file, (setting.output_filename + name_ext).c_str(), "w");
        fprintf(output_file, "Time(s)       Stress(out)        Total_Strain(in)    "
            "   Visco-elastic_Strain    Visco-plastic_Strain \n");

        // Write cable state.
        for (size_t i = 0; i < stressSolver.sigma_cal.size(); i++)
        {
            fprintf(output_file, "% .5E % .5E % .5E % .5E % .5E \n",
                stressSolver.simTime[i], stressSolver.sigma_cal[i], setting.dataIn[i][0],
                stressSolver.eps_ve[i], stressSolver.eps_vp[i]);
        }
        fclose(output_file);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Print log file for error check.
    ////////////////////////////////////////////////////////////////////////////////
    int print_log(string file_name, ErrorCode errCodes, ErrorOut errOut)
    {
        //char str[26];
        //struct tm newtime;
        //std::ofstream log_file;
        //log_file.open(file_name, std::ofstream::app);
        //asctime_s(str, 26, &newtime);
        //log_file << "  " << errOut.message(errCodes) << endl << endl;
        //log_file.close();
        return 0;
    }

} // End of namespace rope.

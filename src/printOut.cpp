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
                strainSolver.simTime[i], setting.dataIn[i][1], strainSolver.eps[i],
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
                stressSolver.simTime[i], stressSolver.sigma_cal[i], setting.dataIn[i][1],
                stressSolver.eps_ve[i], stressSolver.eps_vp[i]);
        }
        fclose(output_file);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Print log file for error check.
    ////////////////////////////////////////////////////////////////////////////////
    int print_log(Setting& setting, ErrorCode errCodes, ErrorOut errOut)
    {
        char buffer[32];
        _time32(&setting.aclock);   // Get time in seconds.
        _localtime32_s(&setting.newtime, &setting.aclock);   // Convert time to struct tm form.

        asctime_s(buffer, 32, &setting.newtime);

        std::ofstream log_file;
        log_file.open(setting.log_filename, std::ofstream::app);
        log_file << "  " << buffer << errOut.message(errCodes) << endl << endl;
        log_file.close();
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Print description to screen.
    ////////////////////////////////////////////////////////////////////////////////
    void print_copyright(void)
    {
        cout << endl << endl;
        cout << "------------------------------------------------------------------"
            "--------------" << endl << endl;
        cout << "           A Nonlinear Synthetic Constitutive Modeling Tool\n";
        fprintf(stdout, " %s Version 1.0     \n", "                         SYNCOM");
        cout << "                  Copyright (c) 2020 Jessica Nguyen\n\n";
        cout << "------------------------------------------------------------------"
            "--------------" << endl << endl;
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// End of Simulation.
    ////////////////////////////////////////////////////////////////////////////////
    void end_simu(void)
    {
        std::ios_base::sync_with_stdio(false);
        std::cout << "\nPress enter to close the application...\n";
        std::cin.ignore();
    }

} // End of namespace rope.

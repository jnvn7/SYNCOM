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


#ifndef SC_stressSolver_api_h
#define SC_stressSolver_api_h

#include "SC_error.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

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

        /// Module selection, time and Stress/Strain input data;
        double tol;
        int limit;
    };

    class stressSolver {

        // Parameters;
        // Nodal properties;
        std::vector<double> te;
        std::vector<double> epsim1;
        std::vector<double> sigmaim1;
        std::vector<double> sigmaim2;
        std::vector<double> g2im1;

        std::vector<double> sigma_yield;
        std::vector<double> eps_vp;
        std::vector<double> sigma_cal;

        std::vector<std::vector<double>> qnim1;

        // Temporary nodal properties;
        std::vector<double> te_Vtemp;
        std::vector<double> sigma_Vtemp;
        std::vector<double> eps_vp_Vtemp;
        std::vector<double> eps_Vtemp;
        std::vector<double> g2_Vtemp;
        std::vector<double> dPsy_Vtemp;

        // Temporary variables;
        int mode, iter;

        double a0, g0, g1, g2, Ep, np, H_vp,
            da0, dg0, dg1, dg2, dEp, dEpm1,  // 1st derivative WRT sigma (applied stress)
            dnp, dnpm1, dH_vp, d2a0, d2g0,   // 2nd derivative WRT sigma (applied stress)
            d2g1, d2g2, d2Ep, d2np, d2H_vp,
            dPsy, d2Psy, d3Psy, Func, DFunc;

        int flag; 

        double jtemp, err, stemp, stemp_new, sumDn1, sumDn2, sumDn3, 
               sumDn4, Exp3, Atemp, Btemp, dAtemp, dBtemp, 
               dCtemp, dExp1, dExp2, dExp3, eps_vp_temp;

        // Functions;
        int calCoeffs(double sigma, double dt);
        int calCoeffs_step(int step_num, double sigma, 
            std::vector<std::vector<double>>& stress_lim, 
            std::vector<double>& xyzCoefs,
            double& xyz, double& dxyz, double& d2xyz);
        void calCoeffs_nostep(std::vector<double>& xyzCoefs, double& xyz,
            double& dxyz, double& d2xyz, double sigma);
        void calDFunc(int mode, int nodeNum, double dt, double te,
                            double sigma, double sigmaim1, double g2im1);
        void calQn(int nodeNum, double sigma, double sigmaim1, double g2, double g2im1, double dPsy);

    public:
        stressSolver(MatProps* mat_props);
        void stressSolver_init(int numNodes);
        ErrorCode validate(void);
        ErrorCode syncom_solver(int nodeNum, double dt, double dataIn, double& stress_SC);
        ErrorCode syncom_init_solver(int nodeNum, double dt, double dataIn, double& stress_SC);
        
        void SC_offInit(void);
        void updateParams(int numNodes, double dt);
        double calStiff(double sigma);
        double get_sigma_yield(int nodeNum) { return sigma_yield[nodeNum]; };
        double get_sigma(int nodeNum) { return sigma_cal[nodeNum]; };
        double get_sigmaim1(int nodeNum) { return sigmaim1[nodeNum]; };
        double get_eps_vp(int nodeNum) { return eps_vp[nodeNum]; };
        double get_eps(int nodeNum) { return eps_Vtemp[nodeNum]; };

        /// Important: The MaterDef.xml should be in the same folder
        std::string sc_inputfile;

        /// The folder is where the MaterDef.xml located
        std::string sc_folder;

        MatProps* material_props;

        /// Path to log file and time parameters;
        std::string log_filename;

#ifndef __unix__
#if _WIN64
        struct tm newtime;
        __time64_t aclock;
#else
        struct tm newtime;
        __time32_t aclock;
#endif
#endif

    };

} // End of namespace rope.

#endif // SC_stressSolver_api_h
#pragma once


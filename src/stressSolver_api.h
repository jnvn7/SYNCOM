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


#ifndef stressSolver_h
#define stressSolver_h

#include "setting.h"
#include "error.h"
#include "printOut_api.h"
#include <math.h>
#include <iostream>

using namespace std;

namespace rope {

    class stressSolver {

        double a0, g0, g1, g2, Ep, np, H_vp,
               da0, dg0, dg1, dg2, dEp, dEpm1,  // 1st derivative WRT sigma (applied stress)
               dnp, dnpm1, dH_vp, d2a0, d2g0,   // 2nd derivative WRT sigma (applied stress)
               d2g1, d2g2, d2Ep, d2np, d2H_vp;  

        double te, dPsy, d2Psy, d3Psy, Func, DFunc,
               epsim1, sigmaim1, sigmaim2, g2im1;

        int mode, iter;

        double sigma_yield, simTime, eps_In, eps_vp, eps_ve, sigma_cal;

        // Temporary variables;
        int flag; 
        double jtemp, err, stemp, stemp_new, sumDn1, sumDn2, sumDn3, 
               sumDn4, Exp3, Atemp, Btemp, dAtemp, dBtemp, 
               dCtemp, dExp1, dExp2, dExp3, eps_vp_temp;

        vector<double> qn;
        vector<double> qnim1;

        int calCoeffs(Setting& setting, double sigma, double dt);
        int calCoeffs_step(int step_num, double sigma, 
            std::vector<std::vector<double>>& stress_lim, 
            std::vector<double>& xyzCoefs,
            double& xyz, double& dxyz, double& d2xyz);
        void calCoeffs_nostep(std::vector<double>& xyzCoefs, double& xyz,
            double& dxyz, double& d2xyz, double sigma);
        void calDFunc(int mode, Setting& setting, double sigma, double dt);
        void calQn(Setting& setting, double sigma);

    public:
        stressSolver(Setting& setting);
        ErrorCode syncom_solver(Setting& setting, double dataIn, double dt);

        double get_sigma_yield(void) { return sigma_yield; };
        double get_simTime(void) { return simTime; };
        double get_sigma(void) { return sigma_cal; };
        double get_eps(void) { return eps_In; };
        double get_eps_ve(void) { return eps_ve; };
        double get_eps_vp(void) { return eps_vp; };

    };

} // End of namespace rope.

#endif // stressSolver_h
#pragma once


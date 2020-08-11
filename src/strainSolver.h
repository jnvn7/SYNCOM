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


#ifndef strainSolver_h
#define strainSolver_h

#include "setting.h"
#include "error.h"
#include <math.h>
#include <iostream>

using namespace std;

namespace rope {

    class strainSolver {

        double a0, g0, g1, g2, Ep, np, H_vp; 
        double te, dPsy, eps_vp_inc, sigmaim1, g2im1;

        double sumDn1, sumDn2, Atemp, Btemp, eps_vp_temp;
        vector<double> qn;
        vector<double> qnim1;

        void calCoeffs(Setting& setting, double sigma, double dt);
        void integrateSR(Setting& setting, double sigma, double dt);
        void calQn(Setting& setting, double sigma);

    public:
        strainSolver(Setting& setting);
        ErrorCode syncom_solver(Setting& setting);

        double sigma_yield;
        vector<double> simTime;
        vector<double> eps;
        vector<double> eps_ve;
        vector<double> eps_vp;
    };

} // End of namespace rope.

#endif // strainSolver_h
#pragma once


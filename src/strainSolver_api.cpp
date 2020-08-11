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

#include "strainSolver_api.h"

namespace rope {
    ///////////////////////////////////////////////////////////////////////////////
    /// Initizlizes Instance;
    ///////////////////////////////////////////////////////////////////////////////
    strainSolver::strainSolver(Setting& setting) {
        
        a0 = g0 = g1 = g2 = Ep = np = H_vp = 0;
        te = dPsy = eps_vp_inc = sigmaim1 = 0;
        g2im1 = 1; sigma_yield = setting.material_props->sigma_yield0;
        
        qn.resize(setting.material_props->lamdaN.size());
        qnim1.resize(setting.material_props->lamdaN.size());
        
        simTime = sigma_In = eps = eps_ve = eps_vp = 0;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Evaluates the material parameters from the input polynomial coefficients;
    /// a0, g0, g1, g2, Ep, np, H;
    ///////////////////////////////////////////////////////////////////////////////
    void strainSolver::calCoeffs(Setting& setting, double sigma) {

        /// Calculate a0, g0, g1, g2, Ep, np, H_vp;
        a0 = g0 = g1 = g2 = Ep = np = H_vp = 0;
        for (size_t i = 0; i < setting.material_props->a0Coefs.size(); i++) {
            a0 += setting.material_props->a0Coefs[i] * pow(sigma, i);
        }

        for (size_t i = 0; i < setting.material_props->g0Coefs.size(); i++) {
            g0 += setting.material_props->g0Coefs[i] * pow(sigma, i);
        }

        for (size_t i = 0; i < setting.material_props->g1Coefs.size(); i++) {
            g1 += setting.material_props->g1Coefs[i] * pow(sigma, i);
        }

        for (size_t i = 0; i < setting.material_props->g2Coefs.size(); i++) {
            g2 += setting.material_props->g2Coefs[i] * pow(sigma, i);
        }

        for (size_t i = 0; i < setting.material_props->EpCoefs.size(); i++) {
            Ep += setting.material_props->EpCoefs[i] * pow(sigma, i);
        }

        for (size_t i = 0; i < setting.material_props->npCoefs.size(); i++) {
            np += setting.material_props->npCoefs[i] * pow(sigma, i);
        }

        for (size_t i = 0; i < setting.material_props->H_vpCoefs.size(); i++) {
            H_vp += setting.material_props->H_vpCoefs[i] * pow(sigma, i);
        }

        /// Calculates dPsy;
        dPsy = 1 / a0 * setting.dt;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Integration of the incremental viscoplastic deformation using quadratic
    ///  Simpson's Rule;
    //////////////////////////////////////////////////////////////////////////////
    void strainSolver::integrateSR(Setting& setting, double sigma) {
        eps_vp_inc = (sigma - setting.material_props->sigma_yield0) / 
                    np * exp(-H_vp / np * te) * setting.dt;
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    /// Calculates previous time step values - Hereditary property Qn
    //////////////////////////////////////////////////////////////////////////////
    void strainSolver::calQn(Setting& setting, double sigma) {

        for (size_t i = 0; i < setting.material_props->lamdaN.size(); i++) {
            qnim1[i] = exp(-setting.material_props->lamdaN[i]*dPsy) * qnim1[i] + 
                        (1 - exp(-setting.material_props->lamdaN[i] * dPsy))/ 
                        (setting.material_props->lamdaN[i] * dPsy) * 
                        (g2 * sigma - g2im1 * sigmaim1);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// SYNCOM - Module 1: Evaluates time history of material behaviors. 
    /// Input: applying stress. Output: time history of strain (deformation).
    //////////////////////////////////////////////////////////////////////////////
    ErrorCode strainSolver::syncom_solver(Setting& setting, double dataIn) {

        /// Calculates instantaneous value for each coefficient;
        calCoeffs(setting, dataIn);
        
        /// Viscoelastic strain;
        sumDn1 = 0; sumDn2 = 0;
        
        for (size_t i1 = 0; i1 < setting.material_props->lamdaN.size(); i1++) {
            sumDn1 = sumDn1 + setting.material_props->Dn[i1] *
                exp(-setting.material_props->lamdaN[i1] * dPsy) * qnim1[i1];
            sumDn2 = sumDn2 + setting.material_props->Dn[i1] *
                (1 - exp(-setting.material_props->lamdaN[i1] * dPsy)) /
                (setting.material_props->lamdaN[i1] * dPsy); 
        }

        Atemp = g0 * setting.material_props->Do +
            g1 * g2 * setting.material_props->sumDn -
            g1 * g2 * sumDn2;
        Btemp = g1 * sumDn1 - g1 * g2im1 * sigmaim1 * sumDn2;

        eps_ve = Atemp * dataIn - Btemp;

        /// Viscoplastic strain;
        if ((dataIn - sigma_yield) > setting.tol && te == 0)
            eps_vp_temp = dataIn / Ep - eps_vp;
        else
            eps_vp_temp = 0;

        // Updates eps_vp_inc;
        if ((dataIn - sigma_yield) > setting.tol) {
            te = te + setting.dt;
            integrateSR(setting, dataIn);
        }
        else {
            eps_vp_inc = 0;
        }

        // Update total viscoplastic deformation;
        eps_vp = eps_vp + eps_vp_inc + eps_vp_temp;

        /// Update total strain;
        eps = eps_ve + eps_vp;

        /// Update sigma_yield, the effective time, and simulation time;
        if ((sigmaim1 - dataIn) > setting.tol && te != 0) {
            if (sigmaim1 > sigma_yield) {
                sigma_yield = sigmaim1;
            }
            te = 0;
        }
        simTime = simTime + setting.dt;
        sigma_In = dataIn;

        // Update previous time step variables;
        calQn(setting, dataIn);  // updates qnim1
        g2im1 = g2;
        sigmaim1 = dataIn;

        if (isnan(eps)) 
            return ErrorCode::NAN_OUTPUT;
        
        return ErrorCode::SUCCESS;

    } // End of syncom_solver

} // End of namespace rope.




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

#include "stressSolver_api.h"

namespace rope {
    ///////////////////////////////////////////////////////////////////////////////
    /// Initizlizes Instance;
    ///////////////////////////////////////////////////////////////////////////////
    stressSolver::stressSolver(Setting& setting) {

        // Initializes zero variables and zero vectors;
        a0 = g0 = g1 = g2 = Ep = np = H_vp = 0;
        da0 = dg0 = dg1 = dg2 = dEp = dnp = dH_vp = 0;
        d2a0 = d2g0 = d2g1 = d2g2 = d2Ep = d2np = d2H_vp = 0;
        te = dPsy = d2Psy = d3Psy = sigmaim1 = sigmaim2 = 0;
        epsim1 = 0;
        g2im1 = 1; sigma_yield = setting.material_props->sigma_yield0;

        qn.resize(setting.material_props->lamdaN.size());
        qnim1.resize(setting.material_props->lamdaN.size());

        simTime = sigma_cal = eps_In = eps_vp = eps_ve = 0;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Evaluates the material parameters from the input polynomial coefficients;
    /// a0, g0, g1, g2, Ep, np, H;
    ///////////////////////////////////////////////////////////////////////////////
    void stressSolver::calCoeffs(Setting& setting, double sigma, double dt) {

        /// Reset values of a0, g0, g1, g2, Ep, np, H_vp, and their derivatives;
        a0 = g0 = g1 = g2 = Ep = np = H_vp = 0;
        da0 = dg0 = dg1 = dg2 = dEp = dnp = dH_vp = 0;
        d2a0 = d2g0 = d2g1 = d2g2 = d2Ep = d2np = d2H_vp = 0;

        /// Calculate a0, g0, g1, g2, Ep, np, H_vp, and their derivatives;
        for (size_t i = 0; i < setting.material_props->a0Coefs.size(); i++) {
            a0 += setting.material_props->a0Coefs[i] * pow(sigma, i);
            if (i > 0) 
                da0 += setting.material_props->a0Coefs[i] * i * pow(sigma, i-1);
            if (i > 1)
                d2a0 += setting.material_props->a0Coefs[i]*i*(i-1)*pow(sigma, i - 2);
        }
       
        for (size_t i = 0; i < setting.material_props->g0Coefs.size(); i++) {
            g0 += setting.material_props->g0Coefs[i] * pow(sigma, i);
            if (i > 0)
                dg0 += setting.material_props->g0Coefs[i] * i * pow(sigma, i - 1);
            if (i > 1)
                d2g0 += setting.material_props->g0Coefs[i]*i*(i - 1)*pow(sigma, i - 2);
        }

        for (size_t i = 0; i < setting.material_props->g1Coefs.size(); i++) {
            g1 += setting.material_props->g1Coefs[i] * pow(sigma, i);
            if (i > 0)
                dg1 += setting.material_props->g1Coefs[i] * i * pow(sigma, i - 1);
            if (i > 1)
                d2g1 += setting.material_props->g1Coefs[i]*i*(i - 1)*pow(sigma, i - 2);
        }

        for (size_t i = 0; i < setting.material_props->g2Coefs.size(); i++) {
            g2 += setting.material_props->g2Coefs[i] * pow(sigma, i);
            if (i > 0)
                dg2 += setting.material_props->g2Coefs[i] * i * pow(sigma, i - 1);
            if (i > 1)
                d2g2 += setting.material_props->g2Coefs[i]*i*(i - 1)*pow(sigma, i - 2);
        }

        for (size_t i = 0; i < setting.material_props->EpCoefs.size(); i++) {
            Ep += setting.material_props->EpCoefs[i] * pow(sigma, i);
            if (i > 0)
                dEp += setting.material_props->EpCoefs[i] * i * pow(sigma, i - 1);
            if (i > 1)
                d2Ep += setting.material_props->EpCoefs[i]*i*(i - 1)*pow(sigma, i - 2);
        }

        for (size_t i = 0; i < setting.material_props->npCoefs.size(); i++) {
            np += setting.material_props->npCoefs[i] * pow(sigma, i);
            if (i > 0)
                dnp += setting.material_props->npCoefs[i] * i * pow(sigma, i - 1);
            if (i > 1)
                d2np += setting.material_props->npCoefs[i]*i*(i - 1)*pow(sigma, i - 2);
        }

        for (size_t i = 0; i < setting.material_props->H_vpCoefs.size(); i++) {
            H_vp += setting.material_props->H_vpCoefs[i] * pow(sigma, i);
            if (i > 0)
                dH_vp += setting.material_props->H_vpCoefs[i] * i * pow(sigma, i - 1);
            if (i > 1)
                d2H_vp += setting.material_props->H_vpCoefs[i]*i*(i - 1)*pow(sigma, i - 2);
        }

        /// Calculates dPsy;
        dPsy = 1 / a0 * dt;
        d2Psy = -pow(a0, -2) * da0 * dt;
        d3Psy = (2*pow(a0, -3)*da0 - pow(a0, -2)*d2a0) * da0 * dt;

        /// Calculates dnpm1;
        dnpm1 = -pow(np, -2) * dnp;

        /// Calculates dEpm1;
        dEpm1 = -pow(Ep, -2) * dEp;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Calculates function's derivatives for Newton-Raphson method;
    ///////////////////////////////////////////////////////////////////////////////
    void stressSolver::calDFunc(int mode, Setting& setting, double sigma, double dt) {

        // Prepares summation terms for construction of Visco-Elastic function;
        sumDn1 = 0; sumDn2 = 0; sumDn3 = 0; sumDn4 = 0;
        for (size_t i = 0; i < setting.material_props->lamdaN.size(); i++) {
            // SumDn1
            sumDn1 += setting.material_props->Dn[i] *
                exp(-setting.material_props->lamdaN[i] * dPsy) *
                qnim1[i];
            
            // SumDn2
            sumDn2 += (setting.material_props->Dn[i] *
                (1 - exp(-setting.material_props->lamdaN[i] * dPsy)) /
                (setting.material_props->lamdaN[i] * dPsy));

            // SumDn3
            dExp1 = -setting.material_props->lamdaN[i] * d2Psy * 
                    exp(-setting.material_props->lamdaN[i] * dPsy);
            sumDn3 += setting.material_props->Dn[i] * dExp1 * qnim1[i];

            // SumDn4
            dExp2 = d2Psy * exp(-setting.material_props->lamdaN[i] * dPsy) / dPsy +
                (1 - exp(-setting.material_props->lamdaN[i] * dPsy)) /
                setting.material_props->lamdaN[i] / dt * da0;
            sumDn4 += setting.material_props->Dn[i] * dExp2;
        }
       
        // Calculates temporary Atemp and Btemp terms;
        Atemp = g0 * setting.material_props->Do + g1 * g2 *
            setting.material_props->sumDn - g1 * g2 * sumDn2;

        Btemp = g1 * sumDn1 - g1 * g2im1 * sigmaim1 * sumDn2;
        
        dAtemp = setting.material_props->Do * dg0 + (dg1 * g2 + g1 * dg2) * 
                    setting.material_props->sumDn - 
                    (dg1 * g2 * sumDn2 + g1 * dg2 * sumDn2 + g1 * g2 * sumDn4);

        dBtemp = g1 * (sumDn3 - g2im1 * sigmaim1 * sumDn4) + 
                    dg1 * (sumDn1 - g2im1 * sigmaim1 * sumDn2);

        // Prepares summation terms for Visco-Plastic function;
        Exp3 = exp(-H_vp / np * te);
        dExp3 = -te * (dH_vp / np + H_vp * dnpm1) * Exp3;

        // Calculates dCd term(Visco-Plastic model);
        if (mode != 0) {
            if (te == dt) {
                dCtemp = 1 / Ep + sigma * dEpm1 + 
                            dt * (1 / np * Exp3 + sigma * 
                                    dnpm1 * Exp3 + sigma * 1 / np * dExp3) - 
                            dt * setting.material_props->sigma_yield0 * 
                            (dnpm1 * Exp3 + 1 / np * dExp3);
            }
            else {
                dCtemp = dt * (1 / np * Exp3 + sigma * 
                            dnpm1 * Exp3 + sigma * 1 / np * dExp3) -
                         dt * setting.material_props->sigma_yield0 * 
                         (dnpm1 * Exp3 + 1 / np * dExp3);
            }
        }
        else {
            dCtemp = 0;
        }
        
        // Calculates DFunc with func = eps - A - B - C;
            DFunc = -Atemp - dAtemp * sigma + dBtemp - dCtemp;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Calculates previous time step values - Hereditary property Qn
    //////////////////////////////////////////////////////////////////////////////
    void stressSolver::calQn(Setting& setting, double sigma) {
        
        for (size_t i = 0; i < setting.material_props->lamdaN.size(); i++) {
            qnim1[i] = exp(-setting.material_props->lamdaN[i] * dPsy) * qnim1[i] +
                (1 - exp(-setting.material_props->lamdaN[i] * dPsy)) /
                (setting.material_props->lamdaN[i] * dPsy) *
                (g2 * sigma - g2im1 * sigmaim1);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// SYNCOM - Module 2: Evaluates time history of material behaviors. 
    /// Input: time history of strain (deformation). Output: applying stress.
    //////////////////////////////////////////////////////////////////////////////
    ErrorCode stressSolver::syncom_solver(Setting& setting, double dataIn, double dt) {

        /////////////////////////////////////////////////////////////////////
        /// VISCO-ELASTIC MODEL ONLY;
        ////////////////////////////////////////////////////////////////////
        // Check if dt is negative;
        if (dt <= 0)
            return ErrorCode::BAD_DT_INPUT;

        if (te == 0 || (sigmaim1 - dataIn) > setting.tol) {
            err = 1; iter = 1; mode = 0;

            // Guesses initial value of stress;
            // Predicts little change in stresses;
            if (abs(sigmaim1 - sigmaim2) < setting.tol)
                stemp = sigmaim1;
            else
                stemp = 2 * sigmaim1 - sigmaim2;

            // Solve for stress iteratively using Visco-Elastic model only;
            while (abs(err) >= setting.tol && iter < setting.limit) {

                /// Calculates instantaneous value for each coefficient;
                calCoeffs(setting, stemp, dt);

                // Updates the function (Func) and its derivative (DFunc);
                calDFunc(mode, setting, stemp, dt);
                Func = dataIn - Atemp * stemp + Btemp - eps_vp;

                // Calculates the new sigma values;
                stemp_new = stemp - Func / DFunc;
                err = stemp_new - stemp;

                if (stemp_new < 0) {
                    stemp = stemp + 0.001;
                    err = 1;
                }
                else if (isnan(stemp))
                    return ErrorCode::NAN_OUTPUT_VISCO_ELASTIC_MODEL;
                else
                    stemp = stemp_new;

                if (abs(stemp_new) < setting.tol) {
                    stemp = 0;
                    err = 0; iter = 1;
                }

                iter = iter + 1;
            } 
        }
        else {
            stemp_new = sigmaim1;
            iter = 5000;

        } // End of Visco-elastic model;

        ////////////////////////////////////////////////////////////////////////////////////
        /// COMBINED VISCO-ELASTIC AND VISCO-PLASTIC MODELS;
        ///////////////////////////////////////////////////////////////////////////////////

        if ((stemp_new - sigma_yield > 1e-6 || iter >= setting.limit)
            && (dataIn - epsim1) >= setting.tol) {

            // Resets conditional variables;
            err = 1; iter = 1; mode = 1;
            te = te + dt;

            // Guesses initial value of stress;
            // Predicts little change in stresses;
            if (abs(sigmaim1 - sigmaim2) < setting.tol)
                stemp = sigmaim1;
            else
                stemp = 2 * sigmaim1 - sigmaim2;
            
            while (abs(err) >= setting.tol && iter <= setting.limit) {
                
                /// Calculates instantaneous value for each coefficient;
                calCoeffs(setting, stemp, dt);

                // Updates visco-plastic strain;
                if (te == dt)
                    eps_vp_temp = stemp / Ep + (stemp - setting.material_props->sigma_yield0) /
                    np * exp(-H_vp / np * te) * dt;
                else
                    eps_vp_temp = eps_vp + (stemp - setting.material_props->sigma_yield0) /
                    np * exp(-H_vp / np * te) * dt;

                // Updates the function (Func) and its derivative (DFunc);
                calDFunc(mode, setting, stemp, dt);
                Func = dataIn - Atemp * stemp + Btemp - eps_vp_temp;

                // Calculates the new sigma values;
                stemp_new = stemp - Func / DFunc;
                err = stemp_new - stemp;

                if (stemp_new < 0) {
                    stemp = stemp + 0.001;
                    iter = 0;
                    err = 1;
                }
                else if (isnan(stemp))
                    return ErrorCode::NAN_OUTPUT_VISCO_ELASTIC_PLASTIC_MODEL;
                else
                    stemp = stemp_new;

                iter = iter + 1;
            }
            
            if (!isnan(stemp)) {
                if (abs(err) <= setting.tol && iter <= setting.limit)
                    eps_vp = eps_vp_temp;
                else
                    return ErrorCode::NO_CONVERGED_SOLUTION;
            }
            else
                return ErrorCode::NAN_OUTPUT_VISCO_ELASTIC_PLASTIC_MODEL; 
        }
        else {
            eps_vp = eps_vp;
        } // End Visco-Elastic and Visco-plastic model;

        /// Updates Current and Previous Time Step Values;
        sigma_cal = stemp_new;
        calQn(setting, sigma_cal);  // updates qnim1
        sigmaim2 = sigmaim1;
        sigmaim1 = sigma_cal;
        g2im1 = g2;

        /// Update sigma_yield, the effective time, and simulation time;
        if ((sigmaim1 - sigma_cal) >
            setting.tol && te > dt) {
            if (sigmaim1 > sigma_yield) {
                sigma_yield = sigmaim1;
            }
            te = 0;
        }

        if ((sigma_yield - sigma_cal) > setting.tol && te == dt) {
            te = 0;
        }

        simTime = simTime + dt;
        eps_ve = dataIn - eps_vp;
        eps_In = dataIn;

        return ErrorCode::SUCCESS;
    }
} // End of namespace rope.




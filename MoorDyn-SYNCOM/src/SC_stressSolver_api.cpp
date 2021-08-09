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

#include "SC_stressSolver_api.h"
#include <iomanip>

namespace rope {
    ///////////////////////////////////////////////////////////////////////////////
    /// Initizlizes Instance and validate input data;
    ///////////////////////////////////////////////////////////////////////////////
    stressSolver::stressSolver(MatProps* mat_props) {
        material_props = mat_props;
        material_props->step_num = std::vector<int>(7, 0);
    };

    ErrorCode stressSolver::validate(void)
    {
        if (material_props->sigma_yield0 < 0
            || material_props->MBL <= 0
            || material_props->Do <= 0) {
            return ErrorCode::BAD_MATERIAL_PROPERTIES_INPUT;
        }
        else {
            return ErrorCode::SUCCESS;
        }

        if (material_props->limit <= 0
            || material_props->tol < 0) {
            return ErrorCode::BAD_NUMERICAL_MATERDEF_INPUT;
        }
        else {
            return ErrorCode::SUCCESS;
        }
    }

    void stressSolver::stressSolver_init(int numNodes) {
        // Initializes zero variables and zero vectors;
        a0 = g0 = g1 = g2 = Ep = np = H_vp = 0;
        da0 = dg0 = dg1 = dg2 = dEp = dnp = dH_vp = 0;
        d2a0 = d2g0 = d2g1 = d2g2 = d2Ep = d2np = d2H_vp = 0;
        dPsy = d2Psy = d3Psy = 0;

        // Nodal properties;
        te.resize(numNodes, 0.0);
        epsim1.resize(numNodes, 0.0);
        sigmaim1.resize(numNodes, 0.0);
        sigmaim2.resize(numNodes, 0.0);
        g2im1.resize(numNodes, 1.0);

        sigma_yield.resize(numNodes, material_props->sigma_yield0);
        eps_vp.resize(numNodes, 0.0);
        sigma_cal.resize(numNodes, 0.0);

        qnim1.resize(numNodes, std::vector<double>(material_props->lamdaN.size(), 0.0));

        te_Vtemp.resize(numNodes, 0.0);
        sigma_Vtemp.resize(numNodes, 0.0);
        eps_vp_Vtemp.resize(numNodes, 0.0);
        eps_Vtemp.resize(numNodes, 0.0);
        g2_Vtemp.resize(numNodes, 1.0);
        dPsy_Vtemp.resize(numNodes, 0.0);
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Evaluates the material parameters from the input polynomial coefficients;
    /// a0, g0, g1, g2, Ep, np, H; With Step Functions;
    ///////////////////////////////////////////////////////////////////////////////
    int stressSolver::calCoeffs_step(int step_num, double sigma,
        std::vector<std::vector<double>>& stress_lim, std::vector<double>& xyzCoefs,
        double& xyz, double& dxyz, double& d2xyz) {

        for (int i = 0; i < step_num; i++) {
            if (i == 0 && sigma <= stress_lim[0][i]) {
                for (int j = 0; j < stress_lim[1][i]; j++) {
                    xyz += xyzCoefs[j] * pow(sigma, j);
                    if (j > 0)
                        dxyz += xyzCoefs[j] * j * pow(sigma, j - 1);
                    if (j > 1)
                        d2xyz += xyzCoefs[j] * j * (j - 1) * pow(sigma, j - 2);
                }
                return 0; // Success
            }
            else if (i == (step_num - 1) && sigma > stress_lim[0][i]) {
                for (int j = stress_lim[1][i]; j < xyzCoefs.size(); j++) {
                    jtemp = j - stress_lim[1][i];
                    xyz += xyzCoefs[j] * pow(sigma, jtemp);
                    if (jtemp > 0)
                        dxyz += xyzCoefs[j] * jtemp * pow(sigma, jtemp - 1);
                    if (jtemp > 1)
                        d2xyz += xyzCoefs[j] * jtemp * (jtemp - 1) * pow(sigma, jtemp - 2);
                }
                return 0; // Success
            }
            else if (sigma > stress_lim[0][i] && sigma <= stress_lim[0][i + 1]) {
                for (int j = stress_lim[1][i]; j < stress_lim[1][i + 1]; j++) {
                    jtemp = j - stress_lim[1][i];
                    xyz += xyzCoefs[j] * pow(sigma, jtemp);
                    if (jtemp > 0)
                        dxyz += xyzCoefs[j] * jtemp * pow(sigma, jtemp - 1);
                    if (jtemp > 1)
                        d2xyz += xyzCoefs[j] * jtemp * (jtemp - 1) * pow(sigma, jtemp - 2);
                }
                return 0; // Success
            }
        }
        return 1; // Failed to find step limits;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Evaluates the material parameters from the input polynomial coefficients;
    /// a0, g0, g1, g2, Ep, np, H; No Step Function;
    ///////////////////////////////////////////////////////////////////////////////
    void stressSolver::calCoeffs_nostep(std::vector<double>& xyzCoefs, double& xyz,
        double& dxyz, double& d2xyz, double sigma) {

        for (int i = 0; i < xyzCoefs.size(); i++) {
            xyz += xyzCoefs[i] * pow(sigma, i);
            if (i > 0)
                dxyz += xyzCoefs[i] * i * pow(sigma, i - 1);
            if (i > 1)
                d2xyz += xyzCoefs[i] * i * (i - 1) * pow(sigma, i - 2);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Evaluates the material parameters from the input polynomial coefficients;
    /// a0, g0, g1, g2, Ep, np, H;
    ///////////////////////////////////////////////////////////////////////////////
    int stressSolver::calCoeffs(double sigma, double dt) {

        /// Reset values of a0, g0, g1, g2, Ep, np, H_vp, and their derivatives;
        a0 = g0 = g1 = g2 = Ep = np = H_vp = 0;
        da0 = dg0 = dg1 = dg2 = dEp = dnp = dH_vp = 0;
        d2a0 = d2g0 = d2g1 = d2g2 = d2Ep = d2np = d2H_vp = flag = 0;
       
        /// Calculate a0, g0, g1, g2, Ep, np, H_vp, and their derivatives;
        // a0;
        if (material_props->step_num[0] == 0) {
            calCoeffs_nostep(material_props->a0Coefs, a0, da0, d2a0, sigma);
        }
        else {
            flag = calCoeffs_step(material_props->step_num[0], sigma, 
                material_props->a0stress_lim,
                material_props->a0Coefs, a0, da0, d2a0);
        }
        if (flag) return flag;

        // g0;
        if (material_props->step_num[1] == 0) {
            calCoeffs_nostep(material_props->g0Coefs, g0, dg0, d2g0, sigma);
        }
        else {
            flag = calCoeffs_step(material_props->step_num[1], sigma, 
                material_props->g0stress_lim,
                material_props->g0Coefs, g0, dg0, d2g0);
        }
        if (flag) return flag;

        // g1;
        if (material_props->step_num[2] == 0) {
            calCoeffs_nostep(material_props->g1Coefs, g1, dg1, d2g1, sigma);
        }
        else {
            flag = calCoeffs_step(material_props->step_num[2], sigma, 
                material_props->g1stress_lim,
                material_props->g1Coefs, g1, dg1, d2g1);
        }
        if (flag) return flag;

        // g2;
        if (material_props->step_num[3] == 0) {
            calCoeffs_nostep(material_props->g2Coefs, g2, dg2, d2g2, sigma);
        }
        else {
            flag = calCoeffs_step(material_props->step_num[3], sigma, 
                material_props->g2stress_lim,
                material_props->g2Coefs, g2, dg2, d2g2);
        }
        if (flag) return flag;

        // Ep;
        if (material_props->step_num[4] == 0) {
            calCoeffs_nostep(material_props->EpCoefs, Ep, dEp, d2Ep, sigma);
        }
        else {
            flag = calCoeffs_step(material_props->step_num[4], sigma, 
                material_props->Epstress_lim,
                material_props->EpCoefs, Ep, dEp, d2Ep);
        }
        if (flag) return flag;

        // np;
        if (material_props->step_num[5] == 0) {
            calCoeffs_nostep(material_props->npCoefs, np, dnp, d2np, sigma);
        }
        else {
            flag = calCoeffs_step(material_props->step_num[5], sigma, 
                material_props->npstress_lim,
                material_props->npCoefs, np, dnp, d2np);
        }
        if (flag) return flag;

        // H_vp;
        if (material_props->step_num[6] == 0) {
            calCoeffs_nostep(material_props->H_vpCoefs, H_vp, dH_vp, d2H_vp, sigma);
        }
        else {
            flag = calCoeffs_step(material_props->step_num[6], sigma, 
                material_props->Hstress_lim,
                material_props->H_vpCoefs, H_vp, dH_vp, d2H_vp);
        }
        if (flag) return flag;

        /// Calculates dPsy;
        dPsy = 1 / a0 * dt;
        d2Psy = -pow(a0, -2) * da0 * dt;
        d3Psy = (2 * pow(a0, -3) * da0 - pow(a0, -2) * d2a0) * da0 * dt;

        /// Calculates dnpm1;
        dnpm1 = -pow(np, -2) * dnp;

        /// Calculates dEpm1;
        dEpm1 = -pow(Ep, -2) * dEp;
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// Calculates function's derivatives for Newton-Raphson method;
    ///////////////////////////////////////////////////////////////////////////////
    void stressSolver::calDFunc(int mode, int nodeNum, double dt, double te, 
                                    double sigma, double sigmaim1, double g2im1) {

        // Prepares summation terms for construction of Visco-Elastic function;
        sumDn1 = sumDn2 = sumDn3 = sumDn4 = Exp3 = dExp3 = 0;
        Atemp = Btemp = dAtemp = dBtemp = dCtemp = DFunc = 0;
        for (size_t i = 0; i < material_props->lamdaN.size(); i++) {
            // SumDn1
            sumDn1 += material_props->Dn[i] *
                exp(-material_props->lamdaN[i] * dPsy) * qnim1[nodeNum][i];
            
            // SumDn2
            sumDn2 += (material_props->Dn[i] *
                (1 - exp(-material_props->lamdaN[i] * dPsy)) /
                (material_props->lamdaN[i] * dPsy));

            // SumDn3
            dExp1 = -material_props->lamdaN[i] * d2Psy * 
                    exp(-material_props->lamdaN[i] * dPsy);
            sumDn3 += material_props->Dn[i] * dExp1 * qnim1[nodeNum][i];

            // SumDn4
            dExp2 = d2Psy * exp(-material_props->lamdaN[i] * dPsy) / dPsy +
                (1 - exp(-material_props->lamdaN[i] * dPsy)) /
                material_props->lamdaN[i] / dt * da0;
            sumDn4 += material_props->Dn[i] * dExp2;
        }
       
        // Calculates temporary Atemp and Btemp terms;
        Atemp = g0 * material_props->Do + g1 * g2 *
            material_props->sumDn - g1 * g2 * sumDn2;

        Btemp = g1 * sumDn1 - g1 * g2im1 * sigmaim1 * sumDn2;
        
        dAtemp = material_props->Do * dg0 + (dg1 * g2 + g1 * dg2) * material_props->sumDn - 
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
                            dt * material_props->sigma_yield0 * 
                            (dnpm1 * Exp3 + 1 / np * dExp3);
            }
            else {
                dCtemp = dt * (1 / np * Exp3 + sigma * 
                            dnpm1 * Exp3 + sigma * 1 / np * dExp3) -
                         dt * material_props->sigma_yield0 * 
                         (dnpm1 * Exp3 + 1 / np * dExp3);
            }
        }
        else {
            dCtemp = 0;
        }
        
        // Calculates DFunc with func = eps - A - B - C;
            DFunc = -Atemp - dAtemp * sigma + dBtemp - dCtemp;

    } // End of calDFunc

    ///////////////////////////////////////////////////////////////////////////////
    /// Calculates previous time step values - Hereditary property Qn
    //////////////////////////////////////////////////////////////////////////////
    void stressSolver::calQn(int nodeNum, double sigma, double sigmaim1, double g2, double g2im1, double dPsy) {
        
        for (size_t i = 0; i < material_props->lamdaN.size(); i++) {
            qnim1[nodeNum][i] = exp(-material_props->lamdaN[i] * dPsy) * qnim1[nodeNum][i] +
                (1 - exp(-material_props->lamdaN[i] * dPsy)) /
                (material_props->lamdaN[i] * dPsy) *
                (g2 * sigma - g2im1 * sigmaim1);
        }
    } // End of calQn

    ///////////////////////////////////////////////////////////////////////////////
    /// Calculates stiffness with no temporal dependence.
    //////////////////////////////////////////////////////////////////////////////
    double stressSolver::calStiff(double sigma) {

        g0 = 0;
        vector<vector<double>> stress_lim = material_props->g0stress_lim;
        if (material_props->step_num[1] == 0) {
            calCoeffs_nostep(material_props->g0Coefs, g0, dg0, d2g0, sigma);
            for (int i = 0; i < material_props->g0Coefs.size(); i++) {
                g0 += material_props->g0Coefs[i] * pow(sigma, i);
            }
        }
        else {
            for (int i = 0; i < material_props->step_num[1]; i++) {
                if (i == 0 && sigma <= stress_lim[0][i]) {
                    for (int j = 0; j < stress_lim[1][i]; j++) {
                        g0 += material_props->g0Coefs[j] * pow(sigma, j);
                    }
                    break; // Success
                }
                else if (i == (material_props->step_num[1] - 1) && sigma > stress_lim[0][i]) {
                    for (int j = stress_lim[1][i]; j < material_props->g0Coefs.size(); j++) {
                        jtemp = j - stress_lim[1][i];
                        g0 += material_props->g0Coefs[j] * pow(sigma, jtemp);
                    }
                    break; // Success
                }
                else if (sigma > stress_lim[0][i] && sigma <= stress_lim[0][i + 1]) {
                    for (int j = stress_lim[1][i]; j < stress_lim[1][i + 1]; j++) {
                        jtemp = j - stress_lim[1][i];
                        g0 += material_props->g0Coefs[j] * pow(sigma, jtemp);
                    }
                    break; // Success
                }
            }
        }       
        double E = 1/(g0 * material_props->Do)*material_props->MBL;
        return E;

    } // End of calStiff

    ///////////////////////////////////////////////////////////////////////////////
    /// SYNCOM_initialize. Use for initialization step in mooring solver
    //////////////////////////////////////////////////////////////////////////////
    ErrorCode stressSolver::syncom_init_solver(int nodeNum, double dt, double dataIn, double& stress_SC) {
     
        /////////////////////////////////////////////////////////////////////
        /// VISCO-ELASTIC MODEL ONLY;
        ////////////////////////////////////////////////////////////////////
        if (dataIn == 0) {
            eps_vp = eps_vp;
            stemp_new = 0;
        }
        else {
            if (te[nodeNum] == 0 || (epsim1[nodeNum] - dataIn) > material_props->tol) {
                err = 1; iter = 1; mode = 0;

                // Guesses initial value of stress;
                // Predicts little change in stresses;
                stemp = sigma_Vtemp[nodeNum];

                // Solve for stress iteratively using Visco-Elastic model only;
                while (abs(err) >= material_props->tol && iter < material_props->limit) {

                    /// Calculates instantaneous value for each coefficient;
                    flag = calCoeffs(stemp, dt);
                    if (flag)
                        return ErrorCode::NON_LOGICAL_COEFFICIENT_INPUT;

                    // Updates the function (Func) and its derivative (DFunc);
                    calDFunc(mode, nodeNum, dt, te[nodeNum], stemp, sigmaim1[nodeNum], g2im1[nodeNum]);
                    Func = dataIn - Atemp * stemp + Btemp - eps_vp[nodeNum];

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

                    if (abs(stemp_new) < material_props->tol) {
                        stemp = 0;
                        err = 0; iter = 1;
                    }

                    iter = iter + 1;
                }
            }
            else {
                stemp_new = sigmaim1[nodeNum];
                iter = material_props->limit;

            } // End of Visco-elastic model;

            if (iter < material_props->limit) {
                stress_SC = stemp_new;
                sigma_Vtemp[nodeNum] = stemp_new;
                eps_Vtemp[nodeNum] = dataIn;
                g2_Vtemp[nodeNum] = g2;
                dPsy_Vtemp[nodeNum] = dPsy;
                return ErrorCode::SUCCESS;
            }
            else 
                return ErrorCode::SYNCOM_INITIALIZATION_FAILED;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// SYNCOM - Module 2: Evaluates time history of material behaviors. 
    /// Input: time history of strain (deformation). Output: applying stress.
    //////////////////////////////////////////////////////////////////////////////
    ErrorCode stressSolver::syncom_solver(int nodeNum, double dt, double dataIn, double& stress_SC) {

        /////////////////////////////////////////////////////////////////////
        /// VISCO-ELASTIC MODEL ONLY;
        ////////////////////////////////////////////////////////////////////
        // Check if dt is negative;
        if (dt <= 0)
            return ErrorCode::BAD_DT_INPUT;
        else if (dataIn < 0)
            return ErrorCode::NEGATIVE_STRAIN_DETECTED;

        if (dataIn <= material_props->tol) {
            stemp_new = 0;
            calCoeffs(stemp_new, dt);
            eps_vp_Vtemp[nodeNum] = eps_vp[nodeNum];
        }
        else {
            if (te[nodeNum] == 0 || (epsim1[nodeNum] - dataIn) > material_props->tol) {
                err = 1; iter = 1; mode = 0;
              
                // Guesses initial value of stress;
                // Predicts little change in stresses;
                if (abs(sigmaim1[nodeNum] - sigmaim2[nodeNum]) < material_props->tol)
                    stemp = sigmaim1[nodeNum];
                else
                    stemp = (2 * sigmaim1[nodeNum] - sigmaim2[nodeNum]);

                // Solve for stress iteratively using Visco-Elastic model only;
                while (abs(err) >= material_props->tol && iter < material_props->limit) {

                    /// Calculates instantaneous value for each coefficient;
                    flag = calCoeffs(stemp, dt);
                    if (flag)
                        return ErrorCode::NON_LOGICAL_COEFFICIENT_INPUT;

                    // Updates the function (Func) and its derivative (DFunc);
                    calDFunc(mode, nodeNum, dt, te[nodeNum], stemp, sigmaim1[nodeNum], g2im1[nodeNum]); 
                    Func = dataIn - Atemp * stemp + Btemp - eps_vp[nodeNum];
                 
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

                    if (abs(stemp_new) < material_props->tol) {
                        stemp = 0;
                        err = 0; iter = 1;
                    }

                    iter = iter + 1; 
                } 
            }
            else {
                stemp_new = sigmaim1[nodeNum];
                iter = material_props->limit;
            } // End of Visco-elastic model;

            ////////////////////////////////////////////////////////////////////////////////////
            /// COMBINED VISCO-ELASTIC AND VISCO-PLASTIC MODELS;
            ///////////////////////////////////////////////////////////////////////////////////

            if ((stemp_new - sigma_yield[nodeNum] > 1e-6 || iter >= material_props->limit)
                && (dataIn - epsim1[nodeNum]) >= material_props->tol) {
                cout << dt << "     " << epsim1[nodeNum] << "     " << dataIn << "     " << DFunc << "    " << "Enter Viscoplastic Solver" << endl;
                // Resets conditional variables;
                err = 1; iter = 1; mode = 1; int iter2 = 1;
                te_Vtemp[nodeNum] = te[nodeNum] + dt;

                // Guesses initial value of stress;
                // Predicts little change in stresses;
                if (abs(sigmaim1[nodeNum] - sigmaim2[nodeNum]) < material_props->tol)
                    stemp = sigmaim1[nodeNum];
                else
                    stemp = 2 * sigmaim1[nodeNum] - sigmaim2[nodeNum];
  
                while (abs(err) >= material_props->tol && iter < material_props->limit) {

                    /// Calculates instantaneous value for each coefficient;
                    flag = calCoeffs(stemp, dt);
                    if (flag)
                        return ErrorCode::NON_LOGICAL_COEFFICIENT_INPUT;

                    // Updates visco-plastic strain;
                    if (te_Vtemp[nodeNum] == dt)
                        eps_vp_temp = stemp / Ep + (stemp - material_props->sigma_yield0) /
                        np * exp(-H_vp / np * te_Vtemp[nodeNum]) * dt;
                    else
                        eps_vp_temp = eps_vp[nodeNum] + (stemp - material_props->sigma_yield0) /
                        np * exp(-H_vp / np * te_Vtemp[nodeNum]) * dt;

                    // Updates the function (Func) and its derivative (DFunc);
                    calDFunc(mode, nodeNum, dt, te_Vtemp[nodeNum], stemp, sigmaim1[nodeNum], g2im1[nodeNum]);
                    Func = dataIn - Atemp * stemp + Btemp - eps_vp_temp;
   
                    // Calculates the new sigma values;
                    stemp_new = stemp - Func / DFunc;
                    err = stemp_new - stemp;

                    if (stemp_new < 0) {
                        stemp = stemp + 0.001;
                        iter = 0; iter2 += 1;
                        err = 1;
                    }
                    else if (isnan(stemp))
                        return ErrorCode::NAN_OUTPUT_VISCO_ELASTIC_PLASTIC_MODEL;
                    else
                        stemp = stemp_new;

                    iter = iter + 1;

                }
         
                if (!isnan(stemp)) {
                    if (abs(err) <= material_props->tol && iter < material_props->limit)
                        eps_vp_Vtemp[nodeNum] = eps_vp_temp;
                    else
                        return ErrorCode::NO_CONVERGED_SOLUTION_VISCO_PLASTIC_SOLVER;
                }
                else
                    return ErrorCode::NAN_OUTPUT_VISCO_ELASTIC_PLASTIC_MODEL;
            }
            else {
                eps_vp_Vtemp[nodeNum] = eps_vp[nodeNum];

            } // End of Visco-Plastic model;

        } // End of Visco-Elastic and Visco-Plastic model;

        // Stores parameters to temporary respository;
        if (iter < material_props->limit) {
            stress_SC = stemp_new;
            sigma_Vtemp[nodeNum] = stemp_new;
            eps_Vtemp[nodeNum] = dataIn;
            g2_Vtemp[nodeNum] = g2;
            dPsy_Vtemp[nodeNum] = dPsy;
        }
        else
            return ErrorCode::NO_CONVERGED_SOLUTION;

        return ErrorCode::SUCCESS;

    } // End of syncom_solver.
     
    ///////////////////////////////////////////////////////////////////////////////
    /// Update time-step variables
    //////////////////////////////////////////////////////////////////////////////
    void stressSolver::updateParams(int numNodes, double dt) {

        for (int nodeNum = 0; nodeNum < numNodes; nodeNum++) {

            /// Transfer data from temporary respository;
            sigma_cal[nodeNum] = sigma_Vtemp[nodeNum];
            te[nodeNum] = te_Vtemp[nodeNum];
            eps_vp[nodeNum] = eps_vp_Vtemp[nodeNum];

            /// Update sigma_cal, sigma_yield, and reset effective time;
            if ((sigmaim1[nodeNum] - sigma_cal[nodeNum]) >
                material_props->tol && te[nodeNum] > dt) {
                if (sigmaim1[nodeNum] > sigma_yield[nodeNum]) {
                    sigma_yield[nodeNum] = sigmaim1[nodeNum];
                }
                te[nodeNum] = 0;
            }

            if ((sigma_yield[nodeNum] - sigma_cal[nodeNum]) > material_props->tol && te[nodeNum] == dt) {
                te[nodeNum] = 0;
            }

            /// Updates Current and Previous Time Step Values;
            calQn(nodeNum, sigma_cal[nodeNum], sigmaim1[nodeNum], g2_Vtemp[nodeNum], g2im1[nodeNum], dPsy_Vtemp[nodeNum]);  // updates qnim1
            sigmaim2[nodeNum] = sigmaim1[nodeNum];
            sigmaim1[nodeNum] = sigma_cal[nodeNum];
            epsim1[nodeNum] = eps_Vtemp[nodeNum];
            g2im1[nodeNum] = g2_Vtemp[nodeNum];
        }

        /// Resets temporary respositories;
        te_Vtemp = te;

    } // End of updateParams

} // End of namespace rope.




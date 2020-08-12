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


#ifndef SynCOM_API_h
#define SynCOM_API_h

#ifndef __unix__
#define DECLDIR __declspec(dllexport)
#else
#define DECLDIR
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    // Initialization.
    int DECLDIR initialize(int module, char input_file[]);
   
    // Main solver - Solve for time history of nonlinear stress/strain development.
    int DECLDIR SynCOM(double dataIn, double dt);

    // Simulation completed;
    int DECLDIR close_app(void);

    // Extracts Total Strain(deformation) result (latest time step ).
    double DECLDIR extract_eps(void);

    // Extracts Visco-elastic Strain(deformation) result (latest time step ).
    double DECLDIR extract_eps_ve(void);

    // Extracts Visco-plastic Strain(deformation) result (latest time step ).
    double DECLDIR extract_eps_vp(void);

    // Extracts Stress result (latest time step ).
    double DECLDIR extract_sigma(void);
/*
    // Clear all global variables and close the program.
    int DECLDIR finish(void); */

#ifdef __cplusplus
}
#endif

#endif // SYNCOM_API_h


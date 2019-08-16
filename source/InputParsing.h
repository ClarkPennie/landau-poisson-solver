/* This is the header file associated to InputParsing.cpp in which the prototypes for the
 * functions contained in that file are declared.  Any other header files which must be linked to
 * for the functions here are also included.
 * 
 */

//************************//
//     INCLUDE GUARDS     //
//************************//

#ifndef INPUTPARSING_H_
#define INPUTPARSING_H_

//************************//
//        INCLUDES        //
//************************//

#include "LP_ompi.h"																					// allows the libraries included, macros defined and external variables declared in LP_ompi.h to be used in the InputParsing functions

//************************//
//   FUNCTION PROTOTYPES  //
//************************//

extern void ReadICOptions(GRVY_Input_Class& iparse);

extern void CheckICOptions(std::string& IC_flag);

extern void ReadICName(GRVY_Input_Class& iparse, std::string IC_flag, std::string& IC_name);

extern void ReadFirstOrSecond(GRVY_Input_Class& iparse);

extern void CheckFirstOrSecond();

extern void ReadGamma(GRVY_Input_Class& iparse, int& gamma);

extern void ReadHomogeneous(GRVY_Input_Class& iparse);

extern void ReadFullandLinear(GRVY_Input_Class& iparse);

extern void ReadLinearLandau(GRVY_Input_Class& iparse);

extern void ReadMassConsOnly(GRVY_Input_Class& iparse);

extern void ReadNoField(GRVY_Input_Class& iparse);

extern void ReadInputParameters(GRVY_Input_Class& iparse, std::string& flag, int& nT,
								int& Nx, int& Nv, int& N, double& nu, double& dt, double& A_amp,
								double& k_wave, double& Lv, double& Lx, double& T_hump, double& shift,
								double& T_0, double& rho_0);

extern void ReadDopingParameters(GRVY_Input_Class& iparse, double& NL, double& NH,
								double& T_L, double& T_R, double& eps);

extern void PrintError(std::string var_name);

#endif /* INPUTPARSING_H_ */

/* This is the source file which contains the subroutines necessary for reading and checking
 * the input data.
 *
 * Functions included: 	ReadICOptions, CheckICOptions, ReadICName, ReadFirstOrSecond,
 * CheckFirstOrSecond, ReadFullandLinear & PrintError.
 *
 *  Created on: Dec 18, 2018
 */

#include "InputParsing.h"																		// InputParsing.h is where the prototypes for the functions contained in this file are declared

void ReadICOptions(GRVY_Input_Class& iparse)													// Function to read the Boolean options for the initital conditions from the input file
{
	// Print a header for the Boolean options from the processor with rank 0:
	if(myrank_mpi==0)																			// only the process with rank 0 will do this
	{
		std::cout << "#=====================================================#" << std::endl;
		std::cout << "#   BOOLEAN OPTIONS TO CHOOSE CURRENT RUN BEHAVIOUR   #" << std::endl;
		std::cout << "#=====================================================#" << std::endl << std::endl;
	}
	
	// Check if each of the IC options have been set and print their values from the 
	// processor with rank 0 (if not, set default value to false):
	if( iparse.Read_Var("Damping",&Damping,false) )													
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> Damping = " << Damping << std::endl;
		}
	}
	if( iparse.Read_Var("TwoStream",&TwoStream,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> TwoStream = " << TwoStream << std::endl;
		}
	}
	if( iparse.Read_Var("FourHump",&FourHump,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> FourHump = " << FourHump << std::endl;
		}
	}
	if( iparse.Read_Var("TwoHump",&TwoHump,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> TwoHump = " << TwoHump << std::endl;
		}
	}
	if( iparse.Read_Var("Doping",&Doping,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> Doping = " << Doping << std::endl;
		}
	}
}

void CheckICOptions(std::string& IC_flag)															// Function to verify the IC options were set correctly
{
	// declare IC_count to count how many ICs have been set
	int IC_count = 0;																				

	// For each IC option, add one to the IC count if it's true and set IC_flag to
	// the relevant location
	if(Damping)
	{
		IC_count++;
		IC_flag.assign("Damping/IC_name");
	}
	if(TwoStream)
	{
		IC_count++;
		IC_flag.assign("TwoStream/IC_name");
	}
	if(FourHump)
	{
		IC_count++;
		IC_flag.assign("FourHump/IC_name");
	}
	if(TwoHump)
	{
		IC_count++;
		IC_flag.assign("TwoHump/IC_name");
	}
	if(Doping)
	{
		IC_count++;
		IC_flag.assign("Doping/IC_name");
	}
	
	// Exit the program if no IC option has been chosen and print a messgage from the
	// processor with rank 0
	if(IC_count == 0)
	{
		if(myrank_mpi==0)
		{
			std::cout << "Program cannot run... No initial condition has been chosen. \n"
					"Please set one of Damping, TwoStream, FourHump or TwoHump to true "
					"in LPsolver-input.txt." << std::endl;
		}
		exit(1);
	}
	// Exit the program if too many IC options are chosen and print a messgage from the
	// processor with rank 0
	if(IC_count > 1)
	{
		if(myrank_mpi==0)
		{
			std::cout << "Program cannot run..." << IC_count << "initial conditions have been chosen." 
						<< std::endl;
			std::cout << "Please ONLY set ONE of Damping, TwoStream, FourHump or TwoHump to true "
							"in LPsolver-input.txt." << std::endl;
		}
		exit(1);
	}
}

void ReadICName(GRVY_Input_Class& iparse, std::string IC_flag, std::string& IC_name)			// Function to read the name of the initital conditions being used from the input file
{
	// Try to read the name of the initial conditions being used for this run
	// (if not available, print a general statement from the process with rank 0)
	if( iparse.Read_Var(IC_flag.c_str(),&IC_name))
	{
		if(myrank_mpi==0)
		{
			std::cout << std::endl << "Using the " << IC_name
					<< " initial conditions for this run." << std::endl << std::endl;
		}
	}
	else
	{
		if(myrank_mpi==0)
		{
			std::cout << std::endl << "The above initial condition which is set equal to 1 "
					"is being used for this run." << std::endl << std::endl;
		}
	}
}

void ReadFirstOrSecond(GRVY_Input_Class& iparse)												// Function to read the Boolean options to decide if this is the first run or not
{
	// Check if each of First or Second have been set and print their values from the
	// processor with rank 0 (if not, set default value to false)
	if( iparse.Read_Var("First",&First,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> First = " << First << std::endl;
		}
	}
	if( iparse.Read_Var("Second",&Second,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> Second = " << Second << std::endl << std::endl;
		}
	}
}

void CheckFirstOrSecond()																		// Function to verify the First or Second options were set correctly
{
	// declare run_count to check if first or second has been set
	int run_count = 0;

	// For each option First or Second, add one to run_count if it's true and print that
	// it's true from the processor with rank 0
	if(First)
	{
		run_count++;
		if(myrank_mpi==0)
		{
			std::cout << "Code is running for the first time with the above initial conditions." 
						<< std::endl << std::endl;
		}
	}
	if(Second)
	{
		run_count++;
		if(myrank_mpi==0)
		{
			std::cout << "Code is picking up from a previous run which should have the above "
					"initial conditions." << std::endl << std::endl;
		}
	}
	// Exit the program more or less than one of these options were set
	if(run_count != 1)
	{
		if(myrank_mpi==0)
		{
			std::cout << "Program cannot run... Need to choose if this is a first run or a subsequent one."
						<< std::endl;
			std::cout << "Please set ONE of First or Second to true in LPsolver-input.txt to decide "
						"if this is the first run or picking up from a previous run." << std::endl;
		}
		exit(1);
	}
}

void ReadFullandLinear(GRVY_Input_Class& iparse)												// Function to read the Boolean option to decide if running with single species collisions or mixed
{
	// Check if FullandLinear has been set and print its value from the 
	// processor with rank 0 (if not, set default value to false):
	if( iparse.Read_Var("FullandLinear",&FullandLinear,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> FullandLinear = " << FullandLinear << std::endl << std::endl;
			if(FullandLinear)
			{
				std::cout << "Running the code with both regular collisions and mixed collisions."
					<< std::endl << std::endl;
			}
			else
			{
				std::cout << "Running with regular single-species collisions." << std::endl << std::endl;
			}
		}
	}
}

void ReadLinearLandau(GRVY_Input_Class& iparse)												// Function to read the Boolean option to decide if running with single species collisions or mixed
{
	// Check if FullandLinear has been set and print its value from the
	// processor with rank 0 (if not, set default value to false):
	if( iparse.Read_Var("LinearLandau",&LinearLandau,false) )
	{
		if(myrank_mpi==0)
		{
			std::cout << "--> LinearLandau = " << LinearLandau << std::endl << std::endl;
			if(LinearLandau)
			{
				std::cout << "Collisions are being modeled by the linear Landau operator (particles collide with those with a Maxwellian density)"
					<< std::endl << std::endl;
			}
			else
			{
				std::cout << "Running with the full Landau collision operator" << std::endl << std::endl;
			}
		}
	}
}

void ReadInputParameters(GRVY_Input_Class& iparse, std::string& flag, int& nT, 
							int& Nx, int& Nv, int& N, double& nu, double& dt, double& A_amp, 
							double& k_wave, double& Lv, double& Lx)								// Function to read all input parameters (IC_flag, nT,  Nx, Nv, N, nu, dt, A_amp, k_wave, L_v & L_x)
{
	// Print a header for the input parameters from the processor with rank 0:
	if(myrank_mpi==0)
	{
		std::cout << "#=====================================================#" << std::endl;
		std::cout << "#           INPUT PARAMETERS FOR CURRENT RUN          #" << std::endl;
		std::cout << "#=====================================================#" << std::endl << std::endl;
	}
	
	// Check if the variable flag has been set and print it from the processor with rank 0
	// (if not, exit)
	if(! iparse.Read_Var("flag",&flag))
	{
		if(myrank_mpi==0)
		{
			std::cout << "Program cannot run..." << std::endl;
			std::cout << "Please set the name of 'flag' in the input file." << std::endl;
		}
		exit(1);
	}
	if(myrank_mpi==0)
	{
		std::cout << "--> Flag associated to this run: " << flag << std::endl << std::endl;
	}

	// Try to read all input parameters to be used for this run and print them from the
	// processor with rank 0 (exit if not available)
	if (! iparse.Read_Var("nT",&nT) )
	{
		PrintError("nT");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %d\n","nT",nT);
	}
	if (! iparse.Read_Var("Nx",&Nx) )
	{
		PrintError("Nx");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %d\n","Nx",Nx);
	}
	if (! iparse.Read_Var("Nv",&Nv) )
	{
		PrintError("Nv");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %d\n","Nv",Nv);
	}
	if (! iparse.Read_Var("N",&N) )
	{
		PrintError("N");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %d\n","N",N);
	}
	if (! iparse.Read_Var("nu",&nu) )
	{
		PrintError("nu");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %g\n","nu",nu);
	}
	if (! iparse.Read_Var("dt",&dt) )
	{
		PrintError("dt");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %g\n","dt",dt);
	}

	// Try to read all input parameters associated to the Damping option 
	// (some have default values but others will exit if not available)
	
	grvy_log_setlevel(GRVY_NOLOG);
	
	if(Damping)
	{
		if (! iparse.Read_Var("Damping/A_amp",&A_amp) )
		{
			PrintError("Damping/A_amp");
			exit(1);
		}
		if (! iparse.Read_Var("Damping/k_wave",&k_wave) )
		{
			PrintError("Damping/k_wave");
			exit(1);
		}
		if (! iparse.Read_Var("Damping/Lv",&Lv) )
		{
			PrintError("Damping/Lv");
			exit(1);
		}
		iparse.Read_Var("Damping/Lx",&Lx,2*PI/k_wave);
	}

	// Try to read all input parameters associated to the TwoStream option 
	// (some have default values but others will exit if not available)
	if(TwoStream)
	{
		if (! iparse.Read_Var("TwoStream/A_amp",&A_amp) )
		{
			PrintError("TwoStream/A_amp");
			exit(1);
		}
		if (! iparse.Read_Var("TwoStream/Lv",&Lv) )
		{
			PrintError("TwoStream/Lv");
			exit(1);
		}
		iparse.Read_Var("TwoStream/Lx",&Lx,2*PI/k_wave);
		k_wave=2*PI/4.;
	}

	// Try to read all input parameters associated to the FourHump option 
	// (some have default values but others will exit if not available)
	if(FourHump)
	{
		iparse.Read_Var("FourHump/A_amp",&A_amp,0.);
		iparse.Read_Var("FourHump/k_wave",&k_wave,0.5);
		if (! iparse.Read_Var("FourHump/Lv",&Lv) )
		{
			PrintError("FourHump/Lv");
			exit(1);
		}
		iparse.Read_Var("FourHump/Lx",&Lx,2*PI/k_wave);
	}

	// Try to read all input parameters associated to the TwoHump option 
	// (some have default values but others will exit if not available)
	if(TwoHump)
	{
		iparse.Read_Var("TwoHump/A_amp",&A_amp,0.);
		iparse.Read_Var("TwoHump/k_wave",&k_wave,0.5);
		if (! iparse.Read_Var("TwoHump/Lv",&Lv) )
		{
			PrintError("TwoHump/Lv");
			exit(1);
		}
		iparse.Read_Var("TwoHump/Lx",&Lx,2*PI/k_wave);
	}

	if(Doping)
	{
		iparse.Read_Var("Doping/A_amp",&A_amp,0.);
		iparse.Read_Var("Doping/k_wave",&k_wave,0.5);
		if (! iparse.Read_Var("Doping/Lv",&Lv) )
		{
			PrintError("Doping/Lv");
			exit(1);
		}
		iparse.Read_Var("Doping/Lx",&Lx,2*PI/k_wave);
	}

	grvy_log_setlevel(GRVY_INFO);

	// Print the values of the parameters read from the specific options from the processor
	// with rank 0
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %g\n","A_amp",A_amp);
		printf("--> %-11s = %g\n","k_wave",k_wave);
		printf("--> %-11s = %g\n","Lv",Lv);
		printf("--> %-11s = %g\n\n","Lx",Lx);
	}
}

void ReadDopingParameters(GRVY_Input_Class& iparse, double& NL, double& NH,
							double& T_L, double& T_R, double& eps)								// Function to read all parameters for a non-uniform doping profile (NL, NH, T_L, T_R, eps)
{
	if(myrank_mpi==0)
	{
		std::cout << "Code is running with a non-uniform background density..." << std::endl;
		std::cout << "Parameters associated to this non-uniform doping profile:" << std::endl << std::endl;
	}

	// Try to read all doping parameters to be used for this run and print them from the
	// processor with rank 0 (exit if not available)
	if (! iparse.Read_Var("Doping/NL",&NL) )
	{
		PrintError("Doping/NL");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %g\n","NL",NL);
	}
	if (! iparse.Read_Var("Doping/NH",&NH) )
	{
		PrintError("Doping/NH");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %g\n","NH",NH);
	}
	if (! iparse.Read_Var("Doping/eps",&eps) )
	{
		PrintError("Doping/eps");
		exit(1);
	}
	if(myrank_mpi==0)
	{
		printf("--> %-11s = %g\n","eps",eps);
	}
	if( iparse.Read_Var("Doping/T_L",&T_L,0.4) )
	{
		if(myrank_mpi==0)
		{
			printf("--> %-11s = %g\n","T_L",T_L);
		}
	}
	if( iparse.Read_Var("Doping/T_R",&T_R,0.4) )
	{
		if(myrank_mpi==0)
		{
			printf("--> %-11s = %g\n","T_R",T_R);
		}
	}
}

void PrintError(std::string var_name)
{
	if(myrank_mpi==0)
	{
		std::cout << "Program cannot run..." << std::endl;
		std::cout << "Please set the value of '" << var_name << "' in the input file." << std::endl;
	}
}


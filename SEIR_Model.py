#!/usr/bin/env python3

""""

    Title :      SEIR_Model
    Notes :      Implementation to evolve the SEIR Model equations developed in Richard Neher's research group at University of Basel. 
    Author:      Shyam Harimohan Menon (Fully original code). Model adapted from Richard Neher's group. See www.covid19-scenarios.org for details. 
    Date  :      11th April 2020

"""

from header import *

# Read user provided arguments
from FunctionInterface import *
ap = argparse.ArgumentParser(description='Command Line Inputs for SEIR_Model.py. All inputs optional. ')
ap.add_argument('-i','--i',metavar='Input_file',default=None,help='Input parameter file for run. By default searches for parameters.ini in current working directory',type=str)
ap.add_argument('-o','--o',metavar='Output_Directory',default=None,help='an argument for specifying the output directory. Default is the current working directory',type=str)
ap.add_argument('-debug','--debug',action='store_true',help='Flag to enable debugging mode, prints extra information to stdout. Default = False. ')

args = vars(ap.parse_args())

if(args['i'] is None) :
	parameter_file = os.getcwd()+ '/parameters.ini'
else : 
	parameter_file= args['i']

if(args['o'] is None) :
	output_directory = os.getcwd()+ '/'
else : 
	output_directory = args['o']


def Equations(t,func,params):


# DESCRIPTION : Returns the derivatives of each of the variables to be solved for in the SEIR Model, to be passed to the ode solver. 

# INPUT:
# 		1. t    : time at which derivative is returned.
# 		2. func : Tuple of the values of the variables at time t. 
# 		3. params : Tuple of the parameters required in the model equations. 

# OUTPUT : 
# 		   Tuple of the derivatives of the variables at time t. 


# Variables to solve for - 
# 1. S : Number of susceptible individuals to the virus
# 2. E : Number of exposed individuals to the virus
# 3. I : Number of infected individuals (when exposed people develop symptoms)
# 4. H : Number of hospitalised individuals 
# 5. C : Number of critical individuals requiring ICU
# 6. R : Nuber of recovered individuals
# 7. D : Number of individuals who died 

	# Values of quantities at time t
	S,E,I,H,C,R,D = np.split(func,7)  
	# Obtain updated value of no_of_age_bins
	from FunctionInterface import no_of_age_bins
	# Reading parameters Start
	#######################################################################################################
	params = np.asarray(params,dtype=np.float64)
	age_dependent_params, age_independent_params = np.split(params,[no_of_age_bins*5])
	t_l , t_i, t_c, t_h, R_0, epsilon, tmax,no_icu_beds,icu_overload_factor = age_independent_params
	N_a,c_a,m_a,f_a,zeta_a = np.split(age_dependent_params,5)
	if(np.sum(C)>no_icu_beds) :
		f_a *= icu_overload_factor 
	# Reading Parameters Over
	#######################################################################################################
	
	# The transmission rate of the virus
	beta = calculate_beta(t,R_0,zeta_a)
	# Account for seasonality if any
	beta *= (1.+(epsilon*np.cos(2*np.pi*(t-tmax))))/t_i 

	#######################################################################################################
	#ODE Equation set

	dS_dt = -1./N_a * beta * S * np.sum(I)
	dE_dt = 1./N_a * beta * S * np.sum(I) - E/t_l
	dI_dt = E/t_l - I/t_i
	dH_dt = (1-m_a)*(I/t_l) + (1-f_a)*C/t_c - H/t_h
	dC_dt = c_a*H/t_h - C/t_c 
	dR_dt = m_a*I/t_i + (1.-c_a)*H/t_h
	dD_dt = f_a*C/t_c

	#######################################################################################################
	derivs = np.concatenate((dS_dt,dE_dt,dI_dt,dH_dt,dC_dt,dR_dt,dD_dt))
	return derivs 

def M(t) : 


# DESCRIPTION : Functional form of the quantitative effect of mitigation measures made by the government. 
# 			  Returned value at time t is the instantaneous factor by which the infection spread rate is decreased. 
# 			  Currently no mitigation implemented, i.e. M(t) = 1.0 for all times.

# INPUT:
# 		1. t    : time at which derivative is returned.
		
# OUTPUT : 
# 		  Instantaneous factor by which the rate is decreased 


	if(t<20.) :
		M = 1.0
	elif((t>20.) and (t<90.)) :
		M = 1.0
	else :
		M = 1.0
	return M


if __name__ == "__main__": 

	
	print("######################################")
	print("Entering COVID 19 Scenario Simulator..")
	print("Reading Parameters....")
	params = set_params()
	
	if(args['debug']) :
		print("\n######################################")
		print("Printing parameters passed by user.")
		print("######################################")
		print_dictionary(params)
	initial_conditions = set_initial_conditions(params['N_a'])
	if(args['debug']) :
		print("######################################")
		print("Printing Initial Conditions")
		print("######################################")
		print_dictionary(initial_conditions)
		print("######################################")

	equation_parameters = list(params['N_a']) + list(params['c_a']) + list(params['m_a']) + list(params['f_a']) + list(params['zeta_a']) \
						+ [params['t_l']] + [params['t_i']] + [params['t_c']] + [params['t_h']]   \
						 + [params['R_0']]  + [params['epsilon']] + [params['tmax']] + [params['no_icu_beds']] + [params['icu_overload_factor']]


	# Read Initial Conditions
	print("Reading Initial Conditions....")
	initial_conditions = get_as_list(initial_conditions)
	initial_conditions = np.asarray(initial_conditions,dtype=np.float64)	

	# Set time axis 
	print("Computing prediction for next {} days...".format(params['end_time']))
	t = np.arange(0,params['end_time'],0.01)

	#######################################################################################################
	# Solve ODE with scipy.integrate.odeint()
	solution = odeint(Equations,initial_conditions.flatten(),t,args=(equation_parameters,),tfirst=True)
	# Convert floats to nearest integer for number of cases
	solution = np.floor(solution)
	# Save solution as a pickle file 
	saveObj(solution,output_directory+'solution')
	print("\nSaved solution as a pickle object in {}".format(output_directory+'solution'))
	Susceptible_Individuals,Exposed_Individuals,Infected_Individuals,Hospitalised_Individuals,Critical_Individuals,Recovered_Individuals,Dead_Individuals = np.split(solution,7,axis=1)
	#Plot the Results
	print("Plotting results...")
	plot_results(Susceptible_Individuals,Exposed_Individuals,Infected_Individuals,Hospitalised_Individuals,
				Critical_Individuals,Recovered_Individuals,Dead_Individuals,params['no_beds'],params['no_icu_beds'],t)
	print("######################################")

	


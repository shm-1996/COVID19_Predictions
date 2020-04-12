#TODO: 
# 1. Read parameters from file.
# 2. Make bins for ages and compute the different fractions for the ages. 
# 3. Account for a multiplicative factor when the threshold for ICU beds is overflowing.
# 4. Fix seasonality Issue
# Assumpttions of the Model -
# 1. Infected + Recovered individuals cannot be affected by the virus again (i.e. herd immunity a possiblity) 
# 2. Only infected people with symptoms (i.e. I subpopulation) can spread the virus, i.e. both asymptomatic and 'in transition' patients cannot spread. 
# 3. The parameters involving fractions of infected populations with severe and mild symptoms (i.e. to hospital or not) are constant in time, which is alright. However the fraction of critical
#.   should be a function of the positive offset from the number of available ICU beds. i.e. f_a = function(C(t)-ICU_beds) which could be linear. They have this in the online version. 


# Questions to ask :
# 1. How long before symptoms show on average?
# 2. How long does the disease last?
# 3. How is the kerala hospital responses, i.e. how severe do cases need to be for an ambulance being called and taken to hospital -> would affect t_i

import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from block_timer.timer import Timer
import matplotlib as mpl
mpl.style.use('classic')
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def Equations(t,func,params):

# Variables to solve for - 
# 1. S : Number of susceptible individuals to the virus
# 2. E : Number of exposed individuals to the virus
# 3. I : Number of infected individuals (when exposed people develop symptoms)
# 4. H : Number of hospitalised individuals 
# 5. C : Number of critical individuals requiring ICU
# 6. R : Nuber of recovered individuals
# 7. D : Number of individuals who died due to COVID-19	

	# Values of quantities at time t
	S,E,I,H,C,R,D = func   
	# Parameters passed
	N, t_l , t_i, t_c, t_h, c_a, m_a, f_a, R_0, zeta, epsilon, tmax = params 
	# The transmission rate of the virus
	beta = calculate_beta(t,R_0,zeta)
	# Account for seasonality if any
	beta *= (1.+(epsilon*np.cos(2*np.pi*(t-tmax))))/t_i 

	dS_dt = -1./N * beta * S * I
	dE_dt = 1./N * beta * S * I - E/t_l
	dI_dt = E/t_l - I/t_i
	dH_dt = (1-m_a)*(I/t_l) + (1-f_a)*C/t_c - H/t_h
	dC_dt = c_a*H/t_h - C/t_c 
	dR_dt = m_a*I/t_i + (1.-c_a)*H/t_h
	dD_dt = f_a*C/t_c

	derivs = [dS_dt,dE_dt,dI_dt,dH_dt,dC_dt,dR_dt,dD_dt]
	return derivs 

#TODO: Read this from a parameter file
def set_params(seasonality) :
	# N : Population of the region
	N = 300000 # 34.8 million people 
	# t_l : Latency time from infection to infectiousness. 
	t_l = 10 # 5 days
	# t_i : time an individual is infectious after which he/she either recovers or falls severely ill, forcing a visit to the hospital.  
	t_i = 3 # 3 days 
	# t_h :  the time a sick person in hospital recovers or deteriorates into a critical state
	t_h = 4 # 4 days 
	# t_c : the time a person remains critical before dying or stabilizing 
	t_c = 5  # 14 days or 2 weeks 
	# m_a : The fraction of infectious that are asymptomatic or mild. 1-m_a corresponds to the fraction that requires hospitalisation. 
	m_a = 0.9 # 90 percent mild cases
	# c_a : the fraction of severe cases (i.e. among hospitalised cases) that turn critical 
	c_a = 0.2 # 20%
	# f_a : the fraction of critical cases that are fatal (i.e. fraction of ICU patients who die)
	f_a = 0.2 # 40%
	# R_0 : Transmittion index. Each individual causes on average R_0 secondary infections while they are infectious.
	R_0 = 2.2 # Riou & Althaus 2020: Pattern of early human-to-human transmission of Wuhan 2019 Novel coronavirus

	# zeta : Degree to which particular age groups are isolated from the rest of the population
	zeta = 1.0

	# no_beds: Number of beds available in the region in hospitals
	no_beds = 5

	#no_icu = Number of Intensive Critical Unit beds available for critical patients (i.e. including ventilators etc)
	no_icu = 1
	# epsilon : Ampilitude of seasonal variation in transmissibility 
	epsilon = 0.8
	# tmax: Time of the year of peak transmission (winter in case of COVID 19)
	tmax = 60.
	
	if(seasonality is False) :
		epsilon = 0.0
	params = N, t_l , t_i, t_c, t_h, c_a, m_a, f_a, R_0, zeta, epsilon, tmax
	return params,no_beds,no_icu
	


# Mitigation measure function. #TODO: Currently no prevention. Assign a functional form to this - they use linear+constant piecewise function 
def M(t) : 
	M = 1.0
	if(t<30) :
		M = 1.0
	elif((t>30) and (t<90)) :
		M = 1.0
	else :
		M = 1.0
	return M

# Return the instantaneous transmission rate of the virus
def calculate_beta(t,R_0,zeta) : 	
	beta = R_0 * zeta * M(t) 
	return beta

#TODO : Read this from a parameter file
def set_initial_conditions(N) :
	# E0: Number of initially exposed individuals (not yet tested positive)
	E0 = 5
	# I0: Number of currently infected individuals (Tested positive)
	I0 = 1
	# H0 : Number of initially hospitalised individuals (due to COVID)
	H0 = 0
	# C0 : Number of initially critical patients of COVID-19 
	C0 = 0
	# Recovered : Number of individuals already recovered from COVID-19
	R0 = 0
	# Number of initially dead patients due to COVID-19 situation 
	D0 = 1

	S0 = N - E0 - I0 - H0-C0 - R0 - D0 

	initial_conditions =  [S0,E0,I0,H0,C0,R0,D0]
	#
	return S0,E0,I0,H0,C0,R0,D0


if __name__ == "__main__": 

	seasonality = False
	plot_directory = "/Users/shm/Desktop/COVID_19/"
	print("######################################")
	print("Entering COVID 19 Scenario Simulator..")
	print("Reading Parameters....")

	
	params, number_of_beds, number_of_icu = set_params(seasonality)
	epsilon, tmax = params[10], params[11]
	if(seasonality) :
		print("Seasonality is switched on with amplitude {} with the peak in {} days".format(epsilon,tmax))
	
	# the latency time for infection
	t_i = params[2]
	#The total population of the region 
	N = params[0] 
	# Initial conditions for the different subsets of the population
	initial_conditions = set_initial_conditions(N)
	
	tstop = 201 # 200 days
	tstep = 1 
	# Set time axis 
	t = np.arange(0,tstop,tstep)
	
	
		
	# Start Timer for Computation Time. 

	#	Solve ODE using scipy
		#solution = scipy.integrate.solve_ivp(fun=Equations,t_span=t_span,y0=initial_conditions,args=(params,beta))	
	solution = odeint(Equations,initial_conditions,t,args=(params,),tfirst=True)
	#solution = solve_ivp(fun=Equations,t_span=(0,tstop),y0=initial_conditions,args=params[:])	
	
	Susceptible_Individuals = solution[:,0]
	Exposed_Individuals = solution[:,1]
	Infected_Individuals = solution[:,2]
	Hospitalised_Individuals = solution[:,3]
	Critical_Individuals = solution[:,4]
	Recovered_Individuals = solution[:,5]
	Dead_Individuals = solution[:,6]	

	plt.clf()
	fig,axs = plt.subplots(ncols = 2,nrows=2,figsize=(12,10))
	axs[0,0].plot(t,Infected_Individuals,color='#7600F590',label='Infected',lw=1.8)
	# axs[0,0].set_title("Infected Individuals",fontsize=20)
	axs[0,1].plot(t,Hospitalised_Individuals,color='#7600F590',label='Hospitalised',lw=1.8)
	axs[0,1].axhline(number_of_beds,color='#FF6D4A90',linestyle='dashed',label='Hospital Beds',lw=1.2)
	axs[1,0].plot(t,Critical_Individuals,color='#7600F590',label='Critical',lw=1.8)
	axs[1,0].axhline(number_of_icu,color='#FF6D4A90',linestyle='dashed',label='ICU Beds',lw=1.2)
	axs[1,1].plot(t,Dead_Individuals,color='#7600F590',label='Dead',lw=1.8)

	axs[0,0].set_ylabel('Infected',fontsize=16)
	axs[0,1].set_ylabel('Hospitalised',fontsize=16)
	axs[1,0].set_ylabel('Critical',fontsize=16)
	axs[1,0].set_xlabel('Time (days)',fontsize=16)
	axs[1,1].set_xlabel('Time (days)',fontsize=16)
	axs[1,1].set_ylabel('Dead',fontsize=16)
	
	shape = np.shape(axs)
	i,j = 0,0
	while i<shape[0] :
		j = 0
		while j<shape[1] :
			# axs[i,j].set_xlabel('Time (days)',fontsize=16)
			# axs[i,j].set_ylabel('Number of People',fontsize=16)
			axs[i,j].xaxis.set_minor_locator(AutoMinorLocator())
			axs[i,j].yaxis.set_minor_locator(AutoMinorLocator())
			axs[i,j].tick_params(axis='x')
			axs[i,j].tick_params(which='major',width=1.5,length=4)
			axs[i,j].tick_params(which='minor',width=0.7,length=2)
			axs[i,j].legend(loc='best')
			axs[i,j].grid(True)
			plt.setp(axs[i,j].spines.values(),linewidth=2.0)
			j +=1
		i += 1	
	
	fig.suptitle("Total Population = {}".format(N),fontsize=24)
	plt.savefig(plot_directory+"results.pdf",bbox_inches = 'tight')


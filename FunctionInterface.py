#!/usr/bin/env python3

""""

    Title :      FunctionInterface
    Notes :      Useful backend functions required for the SEIR_Model.py module. 
    Author:      Shyam Harimohan Menon (Fully original code). Model adapted from Richard Neher's group. See www.covid19-scenarios.org for details. 
    Date  :      11th April 2020

"""


from header import *
import SEIR_Model
no_of_age_bins = 0


def calculate_beta(t,R_0,zeta_a) : 

# DESCRIPTION : Calculates the instantaneous transmission rate beta.
# INPUT:
# 		1. t   : time at which derivative is returned.
# 		2. R_0 : Reproduction no of COVID-19. 
# 		3. zeta_a: Tuple of the isolation factor for age groups.

# OUTPUT : 
# 		   beta : Tuple of the instantaneous transmission rate for all age groups.  

	beta = R_0 * zeta_a * SEIR_Model.M(t) 
	return beta

def set_params() :

# DESCRIPTION : Reads and returns the parameters

# OUTPUT : 
# 		   params: Returns the parameters as a dictionary.  

	# 
	# Read values from each age bin from parameter files
	config = configparser.ConfigParser()
	config.read(SEIR_Model.parameter_file)
	N_population = config['parameters'].getint('no_of_people')

	age_dependent_file = SEIR_Model.output_directory + config['files']['age_dependent_file'] # JSON File

	with open(age_dependent_file) as json_file :
		data = json.load(json_file)

	age_distribution = np.asarray(get_as_list(data['ageDistribution']))
	# SEIR_Model.no_of_age_bins = np.size(age_distribution)
	global no_of_age_bins
	no_of_age_bins = np.size(age_distribution)
	N_a = np.asarray(list(map(int,N_population*age_distribution)))
	m_a = 1. - np.asarray(data['frac']['severe'])
	c_a = np.asarray(data['frac']['critical'])
	f_a = np.asarray(data['frac']['fatal'])
	zeta_a = 1.-np.asarray(data['frac']['isolated'])

	# age independent parameters 
	#Timescaled
	t_l = config['timescales'].getint('symptom_latency_time')
	t_i = config['timescales'].getint('sickness_time')
	t_h = config['timescales'].getint('hospital_time')
	t_c = config['timescales'].getint('critical_time')
	end_time = config['parameters'].getint('end_time')

	#Other important parameters
	R_0 = config['parameters'].getfloat('R_0')
	no_beds = config['parameters'].getint('no_hospital_beds')
	no_icu_beds = config['parameters'].getint('no_icu_beds')
	icu_overload_factor = config['parameters'].getfloat('icu_overload_factor')
	seasonality = config['parameters'].getboolean('seasonality')

	if(seasonality) :
		epsilon = config['seasonality'].getfloat('epsilon')
		tmax    = config['seasonality'].getfloat('tmax')
	else : 
		epsilon = 0.0
		tmax = 0.0

	params = {}
	params = {'N_a':N_a,'c_a':c_a,'m_a':m_a,'f_a':f_a,'zeta_a':zeta_a,'t_l':t_l,
				't_i':t_i,'t_c':t_c,'t_h':t_h,'R_0':R_0,'epsilon':epsilon,
				'tmax':tmax,'no_beds':no_beds,'no_icu_beds':no_icu_beds,'icu_overload_factor':icu_overload_factor,'end_time':end_time}
						
	return params

def set_initial_conditions(N_a) :

# DESCRIPTION : Reads and returns the initial condition.

# OUTPUT : 
# 		   initial_solution: Returns the initial conditions as a dictionary.  

	# Read age dependent parameters from file
	initial_solution = {}
	config = configparser.ConfigParser()
	config.read(SEIR_Model.parameter_file)
	E_0 = config['initial_conditions'].getint('Initial_Exposed')
	I_0 = config['initial_conditions'].getint('Initial_Infected')
	H_0 = config['initial_conditions'].getint('Initial_Hospitalised')
	C_0 = config['initial_conditions'].getint('Initial_Critical')
	R_0 = config['initial_conditions'].getint('Initial_Recovered')
	D_0 = config['initial_conditions'].getint('Initial_Deaths')
	
	# Distribute total cases ~ uniformly across the age bins
	initial_conditions =  [E_0,I_0,H_0,C_0,R_0,D_0]
	# E_0,I_0,H_0,C_0,R_0,D_0 = [distribute_uniformly(i,SEIR_Model.no_of_age_bins) for i in initial_conditions]
	E_0,I_0,H_0,C_0,R_0,D_0 = [distribute_uniformly(i,no_of_age_bins) for i in initial_conditions]
	S_0 = N_a - E_0 - I_0 - H_0 - C_0 - R_0 - D_0 
	initial_solution = {'S_0':S_0,'E_0':E_0,'I_0':I_0,'H_0':H_0,'C_0':C_0,'R_0':R_0,'D_0':D_0}	
	print("\nInitial Total Population : {}".format(np.sum(N_a)))
	print("Initially Infected Individuals : {}".format(np.sum(I_0)))
	print("Initially Exposed Individuals : {}".format(np.sum(E_0)))
	print("Initially Hospitalised Individuals : {}".format(np.sum(H_0)))
	print("Initially Critical Individuals : {}".format(np.sum(C_0)))
	print("Initially Recovered Individuals : {}".format(np.sum(R_0)))	
	print("Initially Dead Individuals : {}".format(np.sum(D_0)))
	print("Initially Susceptible Individuals : {}\n".format(np.sum(S_0)))							
	return initial_solution



def get_as_list(data) :
# DESCRIPTION : Returns all elements of a dictionary as a list.
# INPUT:
#          data: Dictionary to convert to a list.
# OUTPUT : 
# 		   list_data : List obtained from dictionary. 
	list_data = []
	for key in data.keys():
		list_data.append(data[key])	
	return list_data


def distribute_uniformly(value,no_of_bins) :

# DESCRIPTION : Distribute the initial conditions approximately equally to all age groups. First divide value equally across the bins.
#               Then randomly allot the remainder values to random bins.				 
# INPUT:
#          value: Initial condition value.
#          no_of_bins : Number of age bins to distribute the initial value over. 
# OUTPUT : 
# 		   data_array : The initial quantity 
	data_array = np.full((no_of_bins),int(value/no_of_bins))
	
	remainder = value%no_of_bins
	i = 0
	while i<remainder :
		random_no = np.random.randint(0,no_of_bins)
		data_array[random_no] += 1
		i +=1

	return data_array


def print_dictionary(dictionary) :
# DESCRIPTION : Print keys and its values in a dictionary, for keys. 
# INPUT:
#          dictionary : The dictionary for which printing required. 

	for key in dictionary.keys() :
		print("{} = {} \n".format(key,dictionary[key]))

	return


def plot_results(S,E,I,H,C,R,D,number_of_beds,number_of_icu,t) :

# DESCRIPTION : Plot the results in the working directory.  
# INPUT:
#       		S,E,I,H,C,R,D : Solution vector
#				number_of_beds : Number of hospital beds available in the region. 
#				number_of_icu : Number of ICU beds available in the region.
# OUTPUT : 
#				Saves a figure "Cumulative_Cases.pdf" and "New_Cases.pdf" in the working directory. 


	S_total,E_total,I_total,H_total,C_total,R_total,D_total = np.sum(S,axis=1),np.sum(E,axis=1),np.sum(I,axis=1),np.sum(H,axis=1),np.sum(C,axis=1),np.sum(R,axis=1),np.sum(D,axis=1)

	plt.clf()
	fig,axs = plt.subplots(ncols = 2,nrows=2,figsize=(12,10))
	axs[0,0].plot(t,I_total,color='#7600F590',label='Infected',lw=1.8)
	# axs[0,0].set_title("Infected Individuals",fontsize=20)
	axs[0,1].plot(t,H_total,color='#7600F590',label='Hospitalised',lw=1.8)
	axs[0,1].axhline(number_of_beds,color='#FF6D4A90',linestyle='dashed',label='Hospital Beds',lw=1.2)
	axs[1,0].plot(t,C_total,color='#7600F590',label='Critical',lw=1.8)
	axs[1,0].axhline(number_of_icu,color='#FF6D4A90',linestyle='dashed',label='ICU Beds',lw=1.2)
	axs[1,1].plot(t,D_total,color='#7600F590',label='Dead',lw=1.8)

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
			axs[i,j].xaxis.set_minor_locator(AutoMinorLocator())
			axs[i,j].yaxis.set_minor_locator(AutoMinorLocator())
			axs[i,j].tick_params(axis='x')
			axs[i,j].tick_params(which='major',width=1.5,length=4)
			axs[i,j].tick_params(which='minor',width=0.7,length=2)
			axs[i,j].legend(loc='best')
			axs[i,j].grid(True)
			axs[i,j].set_yscale('log')
			plt.setp(axs[i,j].spines.values(),linewidth=2.0)
			j +=1
		i += 1	
	
	#fig.suptitle("Total Population = {}".format(N),fontsize=24)
	plt.tight_layout()
	plt.savefig(SEIR_Model.output_directory+"Total_Cases.pdf",bbox_inches = 'tight')

def saveObj(obj, name):
    """
    Save a pickle object.

    INPUTS:
    ----------
    obj      - the name of the data object to be saved
    name     - the name of the pickle object

    """

    os.system("touch " + name + ".pkl")
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def loadObj(name):
    """
    Load a pickle object.

    INPUTS:
    ----------
    name     - the name of the pickle object

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)



# COVID-19 Predictions
Epidemiological model to predict the impact of COVID-19 on the population of a region. It can predict the number of individuals in a region infected, recovered from the virus, requiring hospitalisation, requiring ICU beds and number of deaths with time. Note that the results of the model are subject to initial conditions and parameters, **and all results should be interpreted with great caution**.

The deterministic compartmental **Susceptible-Exposed-Infected-Recovered (SEIR) Model** was originally developed in [Dr Richard Neher's group](https://neherlab.org/pages/team.html) at the University of Basel. Please do go through their [model](https://covid19-scenarios.org/about) to make sense of and use appropriate parameters for exploring scenarios. This is a rewrite of their model with some mild modifications, primarily to be used for studying the spread and impact of the virus in India.
While the mathematical model has been adapted from their research, **this is fully original code** and no code has been borrowed/adapted from the their [web scenario simulator repository](https://github.com/neherlab/covid19_scenarios).

Note that this is a developer version, i.e. has been created with the intention to modify their base model for specific cases and add some features. Feel free to fork this and make your own changes to the model, however make sure to appropriately credit the original team where relavant.
If however you wish to just explore scenarios, it would be much more effective to use the [interactive simulator in their website](https://covid19-scenarios.org) as it has a user-friendly GUI. However this version could be run if you have python installed in your system.  

To run a scenario with this repository do the following steps-
1. Fork/download this repository from here to your local system.  
2. Convert to an executable by running : `chmod +x main.py`
3. Run the inbuilt scenario as a test by running : `./main.py`
4. Edit parameters.ini to modify the initial conditions or parameters in the model. 
5. Use a different `covid.params.json` for modelling a different state/country. 
6. Run `./main.py -h` for more runtime options. 

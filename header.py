#!/usr/bin/env python



    # Title :      header
    # Notes :      Header file importing packages/modules required for the SEIR model package. 
    # Author:      Shyam Harimohan Menon (Fully original code). Model adapted from Richard Neher's group. See www.covid19-scenarios.org for details. 
    # Date  :      11th April 2020


#Analysis and operations
import numpy as np
from scipy.integrate import odeint

#IO 
import argparse
import configparser
import os
import json
from block_timer.timer import Timer
import pickle

#Plotting
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib import rc
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import ScalarFormatter
mpl.style.use('classic')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',r'\usepackage{mathrsfs}',r'\usepackage{amssymb}']
rc('font', **{'family': 'DejaVu Sans','size':14})
rc('text', usetex=True)





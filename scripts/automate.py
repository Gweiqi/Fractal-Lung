import os
import subprocess
import shutil
import random
import glob
import sys

from multiprocessing import Process
import threading
import time

import numpy as np
from numpy import log, tan, pi, sqrt, exp, sin
from numpy.random import lognormal
import matplotlib.pyplot as plt

from numpy.matlib import repmat
from scipy.signal import cheby2, filtfilt
from scipy.interpolate import BPoly
from scipy.integrate import trapz
from scipy.optimize import minimize

# function to generate input data, flow, tidal volume, time-step, number of breath
def generateFInputSine(number_of_breaths, time_step, breath_period, tidal_volume):

      # number of breaths to be analyzed
      nbr = number_of_breaths

      # time-step (seconds) => inverse of sampling frequency of 'InletFlow'
      dt = time_step

      # breath period (seconds)
      tb = breath_period
      TB = tb*np.ones(nbr)

      # tidal volume (m^3)
      tv = tidal_volume
      TV = tv*np.ones(nbr)

      # slope of flow profile at zero-crossing
      # breath length
      N = np.round(TB/dt)

      # time domain
      time = np.arange(0, nbr*tb, dt)
      if time.size > int(np.round(nbr*tb/dt)):
         time = time[:int(np.round(nbr*tb/dt))]

      #--------------------
      # GENERATE FLOW INPUT
      #--------------------
      # flow profile: simple sine function
      flow = tv*pi/tb*sin(2.*pi/tb*time)

      #--------------------
      # PLOT
      #--------------------
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.plot(time, flow, 'k')
      ax.set_ylabel('flow rate at inlet $Q_\mathrm{in}(t)$ $[m^3/s]$')
      ax.grid('on')
      plt.show()

      #--------------------
      # SAVE
      #--------------------
      cpath = os.getcwd()
      dpath = cpath.replace('scripts','data')

      output = flow
      filename = 'inletFlow'
      np.savetxt(dpath + '/' + filename, output, '%.8f')

      output = np.array([TB, TV, N]).T
      filename = 'TBTVN'
      np.savetxt(dpath + '/' + filename, output, ['%f', '%.8f', '%d'])

      output = np.array([[nbr], [1./dt]]).T
      filename = 'nbfs'
      np.savetxt(dpath + '/' + filename, output, ['%d', '%f'])

# sensitivity analysis function
def sensitivity_function() :
      # plot
      color_index = 0
      fig = plt.figure(); lw = 0.7
      ax_1 = fig.add_subplot(211, ylim=(0,1))
      ax_2 = fig.add_subplot(212)
      ax_1.set_ylabel('concentration $N_2$')
      ax_2.set_ylabel('pleural pressure [Pa]')
      ax_2.set_xlabel('time $t$')

      # define colors for each line plot
      number_of_colors = len(variable_values) + 1
      color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

      # for loop over the range of values
      for x_val in variable_values:

            x_val = float(x_val)
            color_index = color_index + 1

            var = np.ones_like(aw_ind, dtype='f')

            #---------------------
            # change mod. parameter: conctant increase
            #---------------------
            # stretch_width = 0.05
            # phi = -2*stretch_width/N*aw_ind + 1 + stretch_width

            #---------------------
            # change mod. parameter: log-normal distribution
            #---------------------
            # mu, sigma = 0.0, 0.4
            # phi = lognormal(mu, sigma, N)

            #---------------------
            # change mod. parameter: grouped(regional), distributed(random)
            #---------------------
            for k in aw_ind:

                  # either grouped (regional modification) 
                  if (grouped_or_distributed == 'grouped'): 

                        # either compensated pattern         
                        if (compensated_or_noncompensated == 'compensated'):
                            if k < lung_impairment * N :
                              var[k] = x_val
                            if k > (1. - lung_impairment) * N:
                              var[k] = 1./x_val 

                        # or non-compensated pattern  
                        else:
                            if k < lung_impairment * N :
                              var[k] = x_val

                  # or distributed (local modification)
                  else:  

                        # either compensated pattern  
                        if (compensated_or_noncompensated == 'compensated'):
                            if np.mod(k, np.int(1./lung_impairment)) == 0:
                                var[k] *= x_val   
                            elif np.mod(k, np.int(1./lung_impairment)) == 1:
                                var[k] *= 1.0 / x_val  

                        # or non-compensated pattern           
                        else:
                            if np.mod(k, np.int(1./lung_impairment)) == 0:
                                var[k] *= x_val   

            xi = np.ones_like(aw_ind, dtype='f')
            theta = np.ones_like(aw_ind, dtype='f')
            phi = np.ones_like(aw_ind, dtype='f')
            tau = np.ones_like(aw_ind, dtype='f')

            # decide which variable to modify
            if variable_name == 'xi':
                  xi = var
            elif variable_name == 'theta':
                  theta = var
            elif variable_name == 'phi':
                  phi = var
            elif variable_name == 'tau':
                  tau = var
                  
            # generate table for the lung paremter file
            mod_table = np.array([aw_ind, xi, theta, phi, tau]).T

            # save the file as modifyLung.txt, to be read by the exe function
            np.savetxt(data_path + '/modifyLung', mod_table, fmt='%d %f %f %f %f')
            print ('--> Lung Parameter File Modified')
            
            # change the dir to the executive address
            os.chdir(exe_path)
            subprocess.call(['LungModel.exe'])
            print ('--> Processs Finished for : ', variable_name, x_val, grouped_or_distributed, compensated_or_noncompensated)

            # load and unpack data
            print ('--> Loading Data for : ', variable_name, x_val, grouped_or_distributed, compensated_or_noncompensated)
            res_data = np.genfromtxt(data_path + '/primary_results')

            time = res_data[:,0]
            sp1  = res_data[:,1]
            # sp2  = res_data[:,2]
            ppl  = res_data[:,3]

            print ('--> Adding to plot for : ', variable_name, x_val)
            ax_1.plot(time, sp1, color=color[color_index], lw=lw, label = str(variable_name) + ' ' + str(x_val))   
            ax_2.plot(time, ppl, color=color[color_index], lw=lw, label = str(variable_name) + ' ' + str(x_val))
            
            print ('--> Copying Primary Results :', variable_name, x_val, grouped_or_distributed, compensated_or_noncompensated)

            shutil.copyfile(data_path + '/primary_results', data_path 
                                                            + '/primary_results_' 
                                                            + str(variable_name) 
                                                            + '_' 
                                                            + str(x_val) 
                                                            + '_'
                                                            + str(grouped_or_distributed) 
                                                            + '_'
                                                            + str(compensated_or_noncompensated) 
                                                            + '.txt')           
      
      # show the final plot 
      ax_1.legend(loc='best')
      ax_2.legend(loc='best')    
      # plt.show()

      # save
      save_path = data_path.replace('data', 'images')
      plt.savefig(save_path + '/results_' 
                            + str(variable_name) 
                            + str(_) 
                            + str(grouped_or_distributed) 
                            + str(_) 
                            + str(compensated_or_noncompensated) 
                            + '.pdf')

# main function call
# clear the terminal scene
clear = lambda: os.system('cls')
clear()

# current exe and data path
exe_path = os.getcwd().replace('\\scripts','')
data_path = exe_path + '\\data'
print ('exe path is :', exe_path)
print ('data path is :', data_path)

# read input parameters
print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))

nbr = int(sys.argv[1])
print ('number of breath : ', nbr)

grouped_or_distributed = sys.argv[2]
print ('grouped or distributed : ', grouped_or_distributed)

compensated_or_noncompensated = sys.argv[3]
print ('compensated or noncompensated : ', compensated_or_noncompensated)
 
lung_impairment = float(sys.argv[4])
print ('lung_impairment : ', lung_impairment)

variable_name = sys.argv[5]
print ('variable name : ', variable_name)

variable_values = sys.argv[6:]
print ('variable values : ', variable_values)

# length of mod-table. 
# This corresponds with the total number of terminal duct / total number of trumpet lobules, 
# and depends on the value set for d_Lim in the constant/systemProperties file. 
# The total number of terminal ducts is prompted to the console/termial when the flPROG application is run.
N = 304
aw_ind = np.arange(N) # ID of airway unit

#---------------------
# generate flow input file
#---------------------
# time-step (seconds) => inverse of sampling frequency of 'InletFlow'
dt = 0.001
# breath period (seconds)
tb = 3.2
# tidal volume (m^3)
tv = 0.0005

# generate simple sinus flow profile
generateFInputSine(nbr, dt, tb, tv)

# run sensitivity function 
sensitivity_function()

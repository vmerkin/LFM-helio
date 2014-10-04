import os,sys,glob

sys.path.append('../lib')
import util,mas,lfmhlib
import compare_mas_lfm



#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='ConfigScripts/compare.config')
parser.add_argument('PlotConfigFileName',help='The name of the configuration file to use',default='ConfigScripts/plot.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
prm = compare_mas_lfm.params.params(args.ConfigFileName)

#############################
# init variable dictionaires
vars={}
vars['lfm']={}

#############################
### LFM READ STUFF ####
#############################
(vars['lfm']['R'],
vars['lfm']['theta'],
vars['lfm']['phi'],
vars['lfm']['Rc'],
vars['lfm']['thetac'],
vars['lfm']['phic'],
vars['lfm']['br'],
 vars['lfm']['bt'],
 vars['lfm']['bp'],
 vars['lfm']['vr'],
 vars['lfm']['vt'],
 vars['lfm']['vp'],
 vars['lfm']['rho'],
 vars['lfm']['c']) = lfmhlib.read(prm.lfmfile)

vars['lfm']['t'] =  vars['lfm']['c']**2/prm.gamma*1.67e-8/1.38    # FIX ME. UNITS!!!
#############################


#############################
### MAS READ STUFF ####
#############################
vars['mas'] = mas.read_all_vars(prm.masdir,prm.mas_time_label)
mas.set_plot_limits(vars['mas'],args.PlotConfigFileName)
#############################

import compare_mas_lfm.plots

if prm.slice=='k':
    compare_mas_lfm.plots.plot_ij(vars,prm)

else:
    compare_mas_lfm.plots.plot_ik_jk(vars,prm)

    

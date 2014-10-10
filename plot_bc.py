#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('bcFileName',help='BC file to plot')
parser.add_argument('nj',type=int)
parser.add_argument('nk',type=int)
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

#import numpy as np
#import matplotlib.pyplot as plt

from pylab import *
data = loadtxt(args.bcFileName)
"""
figure()
va = sqrt(data[:,2]**2+data[:,3]**2+data[:,4]**2)/sqrt(4.*pi*data[:,8])
vf = sqrt(va**2+data[:,9]**2)
Mf = sqrt(data[:,5]**2+data[:,6]**2+data[:,7]**2)/vf

#pcolormesh(Mf.reshape(args.nj,args.nk))
pcolormesh(data[:,8].reshape(args.nj,args.nk))
colorbar()
show()
"""

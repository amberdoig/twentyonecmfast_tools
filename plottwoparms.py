#this will allow you to pick two items and plot them
#options are all within parm array
#options: z, Nf, Nx, alphaX, Mmin, avg temp

fileglob = '/home/amber/data/alpha0/*'
fileblob = '/home/amber/data/alpha1/*'
#looks at all files within specified alpha folder

from twentyonecmfast_tools import *
#imports everything within directory
 
import matplotlib.pyplot as plt
import numpy as np

#load_andre_models returns 4 sets of data
#this assigns variables to the 4 sets
p,k,d,e=load_andre_models(fileglob)
P,K,D,E=load_andre_models(fileblob)
#p is the parm array
#k is the k array
#d is the delta2 array
#e is the delta2 error array

#p has the option for plotting
#shape of p is 83 data points for each of the 6
#row 0 is z
#row 1 is Nf
#row 2 is Nx
#row 3 is alphaX
#row 4 is Mmin
#row 5 is avg temp

#fixed bottom axis of plot
x=np.argsort(p[:,0])
X=np.argsort(P[:,0])

#fixed data vs other item in same order
plt.plot(p[x,0],p[x,5], color='r', label='alpha0')
plt.plot(P[X,0],P[X,5], color='b', label='alpha1')

#customize plot
plt.ylabel('avg temp')
plt.xlabel('z')
plt.title('z vs. avg temp')
plt.legend(loc='lower right')

#display plot
plt.show()

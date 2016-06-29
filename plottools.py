from twentyonecmfast_tools import *
import matplotlib.pyplot as plt
import numpy as np

def plotter(glob_one,glob_two):
    """User defines two fileglobs
    recommended that variables be set first 
    or whole address can be entered"""

p,k,d,e=load_andre_models(glob_one)
P,K,D,E=load_andre_models(glob_two)

onez=np.argsort(p[:,0])
twoz=np.argsort(P[:,0])

plt.plot(p[onez,0],p[onez,5], color='r', label='glob1')
plt.plot(P[twoz,0],P[twoz,5], color='b', label='glob2')

plt.ylabel('avg temp')
plt.xlabel('z')
plt.title('z vs. avg temp')
plt.legend(loc='lower right')

plt.show()

from twentyonecmfast_tools import *
import matplotlib.pyplot as plt
import numpy as np

def plotter(glob_one,glob_two,y_axis,x_axis=0):
    """User defines two fileglobs
    recommended that variables be set first 
    or whole address can be entered
    0: z
    1: Nf
    2: Nx
    3: alphaX
    4: Mmin
    5: avg temp
    """

    labels_dict = {0 : 'z',
                   1 : 'Nf',
                   2 : 'Nx',
                   3 : 'alphaX',
                   4 : 'Mmin',
                   5 : 'avg temp'}

    p,k,d,e=load_andre_models(glob_one)
    P,K,D,E=load_andre_models(glob_two)

    onez=np.argsort(p[:,x_axis])
    twoz=np.argsort(P[:,x_axis])

    plt.plot(p[onez,x_axis],p[onez,y_axis], color='r', label='glob1')
    plt.plot(P[twoz,x_axis],P[twoz,y_axis], color='b', label='glob2')

    plt.ylabel(labels_dict[y_axis])
    plt.xlabel(labels_dict[x_axis])
    plt.title(labels_dict[y_axis]+' as a function of '+labels_dict[x_axis])
    plt.legend(loc='lower right')

    plt.show()

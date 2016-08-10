"""animation of change in y and z with respect to x
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fileglob='/home/amber/data/binalpha1.5/delta_T_v3_no_halos_z015.00_nf0.999808_useTs1_NX0000.063_alphaX1.5_MminX3.1e+09_aveTb-108.88_Pop2_400_600Mpc'

#set up the figure
fig=plt.figure()  
data=np.fromfile(fileglob,dtype=np.float32)
data=data.reshape((400,400,400))

datamin=np.min(data)
datamax=np.max(data)

xslice=0
im=plt.imshow(data[xslice,:,:],cmap=plt.get_cmap('viridis'),vmin=datamin,vmax=datamax,animated=True)
  
def reset():
    global xslice
    xslice=0
    return im,

def updatefig(*args):
    global data,xslice
    xslice=xslice+1
    im.set_array(data[xslice,:,:])
    return im,

ani=animation.FuncAnimation(fig,updatefig,repeat=True,frames=397,interval=25,blit=True,init_func=reset)
plt.colorbar()
plt.show()

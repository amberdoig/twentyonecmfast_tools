import numpy as np
import matplotlib.pyplot as plt

def sliceplot(file_glob,xslice):
    """User inputs location of binary data file
    and single slice of x axis is returned as plot"""

    data=np.fromfile(file_glob,dtype=np.float32)
    data=data.reshape((400,400,400))
    plt.imshow(data[xslice,:,:])
    plt.colorbar()
    plt.show()    

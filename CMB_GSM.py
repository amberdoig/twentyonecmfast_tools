import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from pygsm import GlobalSkyModel
gsm=GlobalSkyModel()

class CMB_GSM_analysis:
    """description"""
    #init with things that will not change through out.
    #user will have to specify these the first time
    #allow for pass in of cmb file location, but default for ease
    def __init__(self,threshold,sig,freq,nside,LMAX,CMB_file_name='/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits'):
        self.threshold=threshold
        self.sig=sig
        self.freq=freq
        #check that nside value is valid
        if not np.all(hp.isnsideok(nside)):
            raise ValueError("%s is not a valid nside parameter (must be a power of 2, less than 2**30)"%str(nside))
        self.nside=nside
        self.LMAX=LMAX
        self.CMB_file_name=CMB_file_name
        #can read and downgrade cmb in same step
        #use nside instead of 512. change both cmb and gsm to nside
        self.CMB_map=hp.ud_grade(hp.read_map(CMB_file_name),nside)
        #power spectrum of downgraded CMB map
        self.CMB_PS=hp.anafast(self.CMB_map,lmax=LMAX)
        self.ell=np.arange(len(self.CMB_PS))
        self.signal=np.sqrt((np.sum((self.ell*(self.ell+1)*self.CMB_PS)**2))/(len(self.CMB_PS)))
        self.GSM_map=hp.ud_grade(gsm.generate(self.freq),nside)
        self.GSM_PS=hp.anafast(self.GSM_map,lmax=LMAX)

    def make_beam(self,itheta,iphi):
    #calculte the beam based on theta and phi
        #set up grid
        inds=np.arange(12*(self.nside**2))
        #return theta and phi for all other points
        theta,phi=hp.pix2ang(self.nside,inds)
        #return all distances
        angsep=hp.rotator.angdist([itheta,iphi],[theta,phi])
        #return beam
        self.beam=np.exp(-(angsep**2)/(2*self.sig**2))

    #################

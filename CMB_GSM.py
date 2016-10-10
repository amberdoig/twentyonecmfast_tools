import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from pygsm import GlobalSkyModel

def fGSM(freq,scale):
    """fourier of the GSM model 
    used defines freq in MHz such as 5000
    user defines power of 10 scale down such as -4
    """
    gsm=GlobalSkyModel()
    gsmdata=gsm.generate(freq)
    LMAX=1024
    cl=hp.anafast(gsmdata*(10**scale),lmax=LMAX)
    ell=np.arange(len(cl))
    plt.figure()
    scale_str=str(scale)
    freq_str=str(freq)
    plt.plot(ell,ell*(ell+1)*cl,color='b',label='GSM*10^ '+scale_str+' at '+freq_str+' MHz')
    plt.xlabel('ell');plt.ylabel('ell(ell+1)cl(log)');plt.grid()
    hp.write_cl('cl.fits',cl)
    plt.legend(loc='upper right')
    plt.show()

def fCMB():
    """fourier of the CMB data
    """
    wmap_map_I=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    LMAX2=1024
    cl2=hp.anafast(wmap_map_I,lmax=LMAX2)
    ell=np.arange(len(cl2))
    plt.figure() 
    plt.plot(ell,ell*(ell+1)*cl2,color='r',label='CMB') 
    plt.xlabel('ell');plt.ylabel('ell(ell+1)cl');plt.grid()
    hp.write_cl('cl2.fits',cl2)
    plt.legend(loc='upper right')
    plt.show()

def fGandGplusC(freq):
    """GSM and CMB  added together
    compared to GSM alone
    User defines GSM frequency in MHz
    WARNING: Takes some time to run
    Note:
    GSM default shape is 512
    CMB default shape is 2048
    GSM and CMB shape must match to add arrays
    GSM has been cast to a higher nside of 2048
    """
    gsm=GlobalSkyModel()
    gsmdata=gsm.generate(freq)
    gsmdata=hp.ud_grade(gsmdata,2048)
    wmap_map_I=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    LMAX=1024
    cl=hp.anafast(gsmdata+wmap_map_I,lmax=LMAX)
    ell=np.arange(len(cl))
    plt.figure()
    freq_str=str(freq)
    plt.plot(ell,ell*(ell+1)*cl,color='b',label='GSM at '+freq_str+' MHz +CMB')
    plt.xlabel('ell');plt.ylabel('ell(ell+1)cl');plt.grid()
    hp.write_cl('cl.fits',cl)
    gsmdataonly=gsm.generate(freq)
    cl2=hp.anafast(gsmdataonly,lmax=LMAX)
    ell=np.arange(len(cl2))
    plt.plot(ell,ell*(ell+1)*cl2,color='r',label='GSM at '+freq_str+' MHz')
    hp.write_cl('cl2.fits',cl)
    plt.legend(loc='upper right')
    plt.show()

def fGplusC(freq):
    """GSM and CMB  added together
    User defines GSM frequency in MHz
    WARNING:takes some times to run
    Note:
    GSM default shape is 512
    CMB default shape is 2048
    GSM and CMB shape must match to add arrays
    GSM has been cast to a higher nside of 2048
    """
    gsm=GlobalSkyModel()
    gsmdata=gsm.generate(freq)
    gsmdata=hp.ud_grade(gsmdata,2048)
    wmap_map_I=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    LMAX=1024
    cl=hp.anafast(gsmdata+wmap_map_I,lmax=LMAX)
    ell=np.arange(len(cl))
    plt.figure()
    freq_str=str(freq)
    plt.plot(ell,ell*(ell+1)*cl,color='b',label='GSM at '+freq_str+' MHz +CMB')
    plt.xlabel('ell');plt.ylabel('ell(ell+1)cl');plt.grid()
    hp.write_cl('cl.fits',cl)
    plt.legend(loc='upper right')
    plt.show()

def fcrossGandC(freq,scale):
    """GSM+CMB cross with CMB
    compared to CMB only
    User defines GSM freq in MHz 
    User defines scale of GSM such as 0.1 or 0.01
    WARNING: takes some time to run
    """
    gsm=GlobalSkyModel()
    gsmdata=gsm.generate(freq)
    #rather than upgrading GSM, downgrade the CMB to 512
    wmap_map_I=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    CMB=hp.ud_grade(wmap_map_I,512)
    LMAX=1024
    both=(gsmdata*scale)+CMB
#add here: both=both*seebeem
#seebeam=beam()    
    tog=hp.sphtfunc.anafast(both,map2=CMB,lmax=LMAX)
    ell=np.arange(len(tog))
    plt.figure()
    freq_str=str(freq)
    scale_str=str(scale)
    plt.plot(ell,ell*(ell+1)*tog,color='b',label='GSM('+scale_str+') at '+freq_str+' MHZ+CMB cross CMB')
    cl2=hp.anafast(CMB,lmax=LMAX)
    ell=np.arange(len(cl2))
    plt.plot(ell,ell*(ell+1)*cl2,color='r',label='CMB') 
    hp.write_cl('cl2.fits',cl2)
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    hp.write_cl('tog.fits',tog)
    plt.legend(loc='upper right')
    plt.show()

def GandCbeamcross(freq,scale,itheta,iphi,nside,sig=np.pi/10):
    """User defines frequency of GSM and scale to decrease it
    User selects point of interest: theta and phi
    User may designate sigma for beam
    User selects nside
    Default sigma is pi/10
    """
    gsm=GlobalSkyModel()
    gsmdata=gsm.generate(freq)
    #rather than upgrading GSM, downgrade the CMB to 512
    wmap_map_I=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    CMB=hp.ud_grade(wmap_map_I,512)
    LMAX=1024
    both=(gsmdata*scale)+CMB
    showbeam=beam(itheta,iphi,nside,sig)
    npix=12*nside**2
    correction=npix/(np.sum(showbeam))
    both=both*showbeam*correction   
    tog=hp.sphtfunc.anafast(both,map2=CMB,lmax=LMAX)
    ell=np.arange(len(tog))
    plt.figure()
    freq_str=str(freq)
    scale_str=str(scale)
    plt.plot(ell,ell*(ell+1)*tog,color='b',label='GSM('+scale_str+') at '+freq_str+' MHZ+CMB cross CMB')
    cl2=hp.anafast(CMB,lmax=LMAX)
    ell=np.arange(len(cl2))
    plt.plot(ell,ell*(ell+1)*cl2,color='r',label='CMB')
    hp.write_cl('cl2.fits',cl2)
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    hp.write_cl('tog.fits',tog)
    theta_str=str(itheta)
    phi_str=str(iphi)
    plt.legend(loc='upper right')
    plt.title('Point of Interest:('+theta_str+','+phi_str+')')
    plt.ylim(-0.00000001,0.00000006)
    plt.show()


def dist2pt(itheta,iphi,nside):
    """returns array of distances between point of interest and all other points
    """
    #verify that the nside chosen is valid
    check_nside(nside)
    #set up grid
    inds=np.arange(12*(nside**2))
    #return theta and phi for all other points
    theta,phi=hp.pix2ang(nside,inds)
    #return all distances
    return hp.rotator.angdist([itheta,iphi],[theta,phi])

def beam(itheta,iphi,nside,sig=np.pi/10):
    check_nside(nside)
    angsep=dist2pt(itheta,iphi,nside)
    return np.exp(-(angsep**2)/(2*sig**2))

def check_nside(nside):
    """Raises exception if nside is not valid"""
    if not np.all(hp.isnsideok(nside)):
        raise ValueError("%s is not a valid nside parameter (must be a power of 2, less than 2**30)"%str(nside))

"""to plot the showbeam in mollview, remove the comments below and run the file in ipython. you must change the values first"""
#showbeam=beam(np.pi/6,0,512,np.pi/10)
#gsm=GlobalSkyModel()
#gsmdata=gsm.generate(150)
#hp.mollview(showbeam*gsmdata,min=0,max=975)
#plt.show()

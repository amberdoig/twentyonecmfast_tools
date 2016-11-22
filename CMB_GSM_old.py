import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from pygsm import GlobalSkyModel

def GSM_four(freq,scale):
    """plots fourier of the GSM model 
    used defines freq in MHz such as 5000
    user defines power of 10 scale down such as -4
    """
    gsm=GlobalSkyModel()
    #generates the GSM model
    gsmdata=gsm.generate(freq)
    LMAX=1024
    gsmmap=hp.anafast(gsmdata*(10**scale),lmax=LMAX)
    ell=np.arange(len(gsmmap))
    plt.figure()
    #convert variables to text for labels
    scale_str=str(scale)
    freq_str=str(freq)
    plt.plot(ell,ell*(ell+1)*gsmmap,color='b',label='GSM*10^ '+scale_str+' at '+freq_str+' MHz')
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    #delete this line?
    #hp.write_gsmmap('gsmmap.fits',gsmmap)
    plt.legend(loc='upper right')
    plt.show()

def CMB_four():
    """plots fourier of the CMB data
    """
    #the raw CMB data
    cmbdata=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    LMAX2=1024
    cmbmap=hp.anafast(cmbdata,lmax=LMAX2)
    ell=np.arange(len(cmbmap))
    plt.figure() 
    plt.plot(ell,ell*(ell+1)*cmbmap,color='r',label='CMB') 
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    #delete this line?
    #hp.write_cmbmap('cmbmap.fits',cmbmap)
    plt.legend(loc='upper right')
    plt.show()

def G_plus_C_vs_G_four(freq):
    """plots fourier of GSM and CMB  added together
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
    #generates a model of the GSM
    gsmdata=gsm.generate(freq)
    #upgrade the GSM data to 2048 
    gsmdata=hp.ud_grade(gsmdata,2048)
    #the raw CMB data
    cmbdata=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    LMAX=1024
    bothmap=hp.anafast(gsmdata+cmbdata,lmax=LMAX)
    ell=np.arange(len(bothmap))
    plt.figure()
    freq_str=str(freq)
    plt.plot(ell,ell*(ell+1)*bothmap,color='b',label='GSM at '+freq_str+' MHz +CMB')
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    #delete this line?
    #hp.write_bothmap('bothmap.fits',bothmap)
    gsmdataonly=gsm.generate(freq)
    gsmmaponly=hp.anafast(gsmdataonly,lmax=LMAX)
    ell=np.arange(len(gsmmaponly))
    plt.plot(ell,ell*(ell+1)*gsmmaponly,color='r',label='GSM at '+freq_str+' MHz')
    #delete this line?
    #hp.write_gsmmaponly('gsmmaponly.fits',gsmmaponly)
    plt.legend(loc='upper right')
    plt.show()

def G_plus_C_four(freq):
    """plots fourier of GSM and CMB  added together
    User defines GSM frequency in MHz
    WARNING:takes some times to run
    Note:
    GSM default shape is 512
    CMB default shape is 2048
    GSM and CMB shape must match to add arrays
    GSM has been cast to a higher nside of 2048
    """
    gsm=GlobalSkyModel()
    #creates a model of the GSM
    gsmdata=gsm.generate(freq)
    #upgrades the GSM to 2048
    gsmdata=hp.ud_grade(gsmdata,2048)
    cmbdata=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    LMAX=1024
    bothmap=hp.anafast(gsmdata+cmbdata,lmax=LMAX)
    ell=np.arange(len(bothmap))
    plt.figure()
    freq_str=str(freq)
    plt.plot(ell,ell*(ell+1)*bothmap,color='b',label='GSM at '+freq_str+' MHz +CMB')
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    #delete this line?
    #hp.write_bothmap('bothmap.fits',bothmap)
    plt.legend(loc='upper right')
    plt.show()

def G_plus_C_cross_C_four(freq,scale):
    """plot of fourier of
    GSM(scale)+CMB crossed with CMB
    compared to CMB only
    User defines GSM freq in MHz 
    User defines scale of GSM such as 0.1 or 0.01
    WARNING: takes some time to run
    """
    gsm=GlobalSkyModel()
    #create the GSM model
    gsmdata=gsm.generate(freq)
    #rather than upgrading GSM, downgrade the CMB to 512
    cmbdata=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    cmbmap=hp.ud_grade(cmbdata,512)
    LMAX=1024
    both=(gsmdata*scale)+cmbmap 
    #cross cmb+gsm with the cmb
    c_plus_g_cross_c=hp.sphtfunc.anafast(both,map2=cmbmap,lmax=LMAX)
    ell=np.arange(len(c_plus_g_cross_c))
    plt.figure()
    #format variables as text for labels
    freq_str=str(freq)
    scale_str=str(scale)
    plt.plot(ell,ell*(ell+1)*c_plus_g_cross_c,color='b',label='GSM('+scale_str+') at '+freq_str+' MHZ+CMB cross CMB')
    cmbmaponly=hp.anafast(cmbdata,lmax=LMAX)
    ell=np.arange(len(cmbmaponly))
    plt.plot(ell,ell*(ell+1)*cmbmaponly,color='r',label='CMB') 
    #delete this line?
    #hp.write_cmbmaponly('cmbmaponly.fits',cmbmaponly)
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    #delete this line?
    #hp.write_cp_plus_g_cross_c('c_plus_g_cross_c.fits',c_plus_g_cross_c)
    plt.legend(loc='upper right')
    plt.show()

def Beam_G_plus_C_cross_C_four(freq,scale,itheta,iphi,nside,sig=np.pi/10):
    """plot of fourier for
    GSM(scale)+CMB crossed with CMB
    beam applied then accounted for
    all compared to CMB alone
    User defines frequency of GSM and scale to decrease it
    User selects point of interest: theta and phi
    User may designate sigma for beam
    User selects nside
    Default sigma is pi/10
    """
    gsm=GlobalSkyModel()
    #generates GSM model
    gsmdata=gsm.generate(freq)
    #rather than upgrading GSM, downgrade the CMB to 512
    cmbdata=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    cmbmap=hp.ud_grade(cmbdata,512)
    LMAX=1024
    both=(gsmdata*scale)+cmbmap
    #run the beam function with the given input
    showbeam=beam(itheta,iphi,nside,sig)
    npix=12*nside**2
    #correct the data to account for the beam
    correction=npix/(np.sum(showbeam))
    both=both*showbeam*correction   
    #cross GSM+CMB with CMB and transform
    c_plus_g_cross_c=hp.sphtfunc.anafast(both,map2=cmbmap,lmax=LMAX)
    ell=np.arange(len(c_plus_g_cross_c))
    plt.figure()
    #convert variables into text for labels
    freq_str=str(freq)
    scale_str=str(scale)
    plt.plot(ell,ell*(ell+1)*c_plus_g_cross_c,color='b',label='GSM('+scale_str+') at '+freq_str+' MHZ+CMB cross CMB')
    #transform CMB
    cmbmaponly=hp.anafast(cmbmap,lmax=LMAX)
    ell=np.arange(len(cmbmaponly))
    plt.plot(ell,ell*(ell+1)*cmbmaponly,color='r',label='CMB')
    #delete this line?
    #hp.write_cmbmaponly('cmbmaponly.fits',cmbmaponly)
    plt.xlabel('ell');plt.ylabel('ell(ell+1)');plt.grid()
    #delete this line?
    #hp.write_c_plus_g_cross_c('c_plus_g_cross_c.fits',c_plus_g_cross_c)
    #convert variable into text for labels
    theta_str=str(itheta)
    phi_str=str(iphi)
    plt.legend(loc='upper right')
    plt.title('Point of Interest:('+theta_str+','+phi_str+')')
    #set the limits for the y axis
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



def G_plus_C_cross_C_beam_chi(freq,scale,itheta,iphi,nside,sig,cmbmap,cmbmaponly):
    """returns the ratio of chi/signal
    for beam applied to GSM(scale)+CMB 
    crossed with CMB
    versus CMB alone
    uses same input as Beam_G_plus_C_cross_C_four()
    but returns value of ratio only, no plot
    CANNOT BE RUN ALONE, REFERENCES ANOTHER FUNCTION
    """
    gsm=GlobalSkyModel()
    #generates GSM model
    gsmdata=gsm.generate(freq)
    #rather than upgrading GSM, downgrade the CMB to 512
    
    #following lines moved into find chi function
    #cmbdata=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    #cmbmap=hp.ud_grade(cmbdata,512)
    
    LMAX=1024
    both=(gsmdata*scale)+cmbmap
    #run the beam function with the given input
    showbeam=beam(itheta,iphi,nside,sig)
    npix=12*nside**2
    #correct the data to account for the beam
    correction=npix/(np.sum(showbeam))
    both=both*showbeam*correction
    #cross GSM+CMB with CMB and transform
    c_plus_g_cross_c=hp.sphtfunc.anafast(both,map2=cmbmap,lmax=LMAX)
    
    #moved to find function
    #cmbmaponly=hp.anafast(cmbmap,lmax=LMAX)
    
    ell=np.arange(len(cmbmaponly)) 
    #find chi squared, the deviation from the desired data
    difference=c_plus_g_cross_c-cmbmaponly
    chi_squared=(np.sum((ell*(ell+1)*difference)**2))/(len(c_plus_g_cross_c))
    signal_squared=(np.sum((ell*(ell+1)*cmbmaponly)**2))/(len(c_plus_g_cross_c))
    chi=np.sqrt(chi_squared)
    signal=np.sqrt(signal_squared)
    ratio=chi/signal
    return ratio


def find_scale_G_plus_C_cross_C_beam(threshold,freq,itheta,iphi,nside,sig=np.pi/10):
    """runs through G_plus_C_cross_C_beam_chi() with decreasing scaled
    until ration is equal to or less than threshold
    """
    #ratio must start larger than we want
    ratio=threshold*2
    #initial scale value must be set
    #must account for first decrease
    scale=0.0009765625
    cmbdata=hp.read_map('/home/amber/data/cmb/COM_CMB_IQU-commander-field-Int_2048_R2.01_full.fits')
    cmbmap=hp.ud_grade(cmbdata,512)
    LMAX=1024
    cmbmaponly=hp.anafast(cmbmap,lmax=LMAX)
    while (ratio > threshold):
        #decrease the existing scale
        scale=scale*.25
        #run the function again with the new scale and return the ratio
        ratio=G_plus_C_cross_C_beam_chi(freq,scale,itheta,iphi,nside,sig,cmbmap,cmbmaponly)
        #ratio will retur to the beginning of while loop. 
        #if larger than threshold will run again
        #if same or less, will return the current scale value
    return scale




"""to plot the showbeam in mollview, remove the comments below and run the file in ipython. you must change the values first"""
#showbeam=beam(np.pi/6,0,512,np.pi/10)
#gsm=GlobalSkyModel()
#gsmdata=gsm.generate(150)
#hp.mollview(showbeam*gsmdata,min=0,max=975)
#plt.show()


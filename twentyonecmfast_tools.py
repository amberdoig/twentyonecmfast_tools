import numpy as n,os
from glob import glob
from scipy.interpolate import LinearNDInterpolator,interp1d
from capo import cosmo_units
from scipy import integrate
import matplotlib.pyplot as plt

def build_model_interp(parm_array,delta2_array,k_array,redshift,
    regrid_ks=None):
    #input an array of models, the parameters for each model, list of k modes
    #   and the desired redshift
    #return a list of ks and a list of matching interpolation functions,
    # function has call f(log10(Nx),alphaX,log10(MminX))
    #parm_array expected to be nmodels,nparms
    #with columns (z,Nf,Nx,alphaX,MminX,other-stuff....)
    #delta2_array expected to be nmodels,nkmodes
    #NOTE: assumes all models are computed at the same k modes
    closest_redshift = parm_array[n.abs(parm_array[:,0]-redshift).argmin(),0]
    model_points = parm_array[parm_array[:,0]==closest_redshift,2:5]

    #interpolate NX and Mmin in log space
    model_points[:,0] = n.log10(model_points[:,0])
    model_points[:,2] = n.log10(model_points[:,2])
    #get the power spectrum values that go with this redshift
    raw_model_values = delta2_array[parm_array[:,0]==closest_redshift]
    model_values = [[]]*len(raw_model_values)
    if not regrid_ks is None: #regrid to a different set of k bins
        for i in xrange(len(raw_model_values)):
            model_values[i] = interp1d(k_array,raw_model_values[i,:])(regrid_ks)
        model_values = n.array(model_values)
        k_array = regrid_ks
        print "interpolated sim shape after regridding",model_values.shape
    else:
        model_values = raw_model_values.copy()
    #for a single redshift, build an interplation for each k mode
    Pk_models_atz = []
    for ki,k in enumerate(k_array):

        M = LinearNDInterpolator(model_points,model_values[:,ki])
        Pk_models_atz.append(M)
    return Pk_models_atz
def build_tau_interp_model(parm_array):
    #interpolate NX and Mmin in log space
    alphaXs = n.sort(list(set(parm_array[:,3])))
    Mmins = n.sort(list(set(parm_array[:,4])))
    Nxs = n.sort(list(set(parm_array[:,2])))
    taus = []
    for Nx in Nxs:
        for alphaX in alphaXs:
            for Mmin in Mmins:
                _slice = n.argwhere(all_and([
                                    parm_array[:,2]==Nx,
                                    parm_array[:,3]==alphaX,
                                    parm_array[:,4]==Mmin]
                                    ))
                taus.append([n.log10(Nx),alphaX,n.log10(Mmin),
                            nf_to_tau(parm_array[_slice,0].squeeze(),
                            parm_array[_slice,1].squeeze())])
    taus = n.array(taus)
    return LinearNDInterpolator(taus[:,:3],taus[:,3])

def all_and(arrays):
    #input a list or arrays
    #output the arrays anded together
    if len(arrays)==1:return arrays
    out = arrays[0]
    for arr in arrays[1:]:
        out = n.logical_and(out,arr)
    return out

def load_andre_models(fileglob):
    #input a string that globs to the list of input model files
    #return arrays of parameters,k modes, delta2,and delt2 error
    #parm_array expected to be nmodels,nparms
    #with columns (z,Nf,Nx,alphaX,Mmin,other-stuff....)
    #delta2_array expected to be nmodels,nkmodes
    filenames = glob(fileglob)
    parm_array = []
    k_array = []
    delta2_array = []
    delta2_err_array = []
    for filename in filenames:
        parms = os.path.basename(filename).split('_')
        if parms[0].startswith('reion'):continue
        parm_array.append(map(float,[parms[3][1:],
                            parms[4][2:], #Nf
                            parms[6][2:], #Nx
                            parms[7][-3:], #alphaX
                            parms[8][5:], #Mmin
                            parms[9][5:]]))
        D = n.loadtxt(filename)
        k_array.append(D[:,0])
        delta2_array.append(D[:,1])
        delta2_err_array.append(D[:,2])
    parm_array = n.array(parm_array)
    raw_parm_array = parm_array.copy()
    k_array = n.ma.array(k_array)
    raw_k_array = k_array.copy()
    delta2_array = n.ma.masked_invalid(delta2_array)
    raw_delta2_array = delta2_array.copy()
    delta2_err_array = n.ma.array(delta2_err_array)
    return parm_array,k_array,delta2_array,delta2_err_array
def load_andre_global_models(fileglob):
    #input a string that globs to the list of input model files
    #return a concatenated array of parameters
    #columns Nx,alphaX,Mmin
    # and global histories
    #dimensions len(parms) x n_redshifts x nparms
    #there are other columns that I don't know what they are
    filenames = glob(fileglob)
    parm_array = []
    global_evolution = []
    for filename in filenames:
        if not os.path.basename(filename).startswith('global'): continue
        parms = os.path.basename(filename).split('_')
        parm_array.append(map(float,
                            [parms[5][2:], #Nx
                             parms[6][-3:], #alphaX
                             parms[7][5:], #Mmin
                            ]))
        D = n.loadtxt(filename)
        global_evolution.append(D)
    return n.array(parm_array),global_evolution

def nf_to_tau(z,nf):
    #based on Liu et al 1509.08463
    """
    i) Take your ionization history, x_{HII} (z).
    ii) Numerically compute the integral \int_0^{zCMB} dz x_{HII} (1+z)^2 / sqrt{OmL + Omm (1+z)^3}.
    I would take OmL = 1 - Omm and Omm = 0.3089.
    In practice the integral doesn't need to be literally taken to zCMB since the ionization fraction is basically zero well before you hit those redshifts
    iii) Multiply by 0.00210228. This includes all the constants in Eq. (10) of my paper.
    iv) Add 0.001. That accounts for helium reionization (if you want to be more precise, it's 0.001223 for helium reionization happening at z = 3.5 and 0.000986 for z = 3.0).
    v) You have tau_CMB!
    """
    Omm = 0.3089
    Oml = 1 - Omm
    coeff = 0.00210228 #this includes all the constants out front of eq 10
    z = n.concatenate([n.linspace(0,z.min()),z])
    nf = n.concatenate([n.zeros(50),nf])
    xHI = interp1d(z,1-nf)
    E = lambda z: xHI(z) * (1+z)**2 / n.sqrt(Oml + Omm * (1+z)**3)
    tau_H  = integrate.quad(E,z.min(),z.max())[0]*coeff
    tau_He = 0.001223 #That accounts for helium reionization
    tau = tau_H + tau_He
    return tau


def compare_runs(*args):
    # Input list of file globs to 21cmfast runs
    # Will generate some plots to compare the runs
    nruns = len(args)
    parms, ks, delta2s, errs = [], [], [], []
    for run in args:
        temp = load_andre_models(run)
        parms.append(temp[0])
        ks.append(temp[1])
        delta2s.append(temp[2])
        errs.append(temp[3])

    plt.figure('comparison')
    for run in xrange(nruns):
        plot(parms[run][:, 0], parms[run][:, 5])
        

import numpy as np
from mpi4py import MPI
import time
import pathlib
from dedalus import public as de

import logging
logger = logging.getLogger(__name__)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def spec1d_classic(domain,fields):
    """
    [modn,En] = spec1d_classic(domain,[u,v,w]) or
    [modn,En] = spec1d_classic(domain,[u,v])

    input: domain, fields (list of dedalus fields, like [u,v,w])
    output: [modn,En], 1-dimensional power spectrum of u['c']**2+v['c']**2 + ...
    modn is the mode number magnitude.
    
    Notes: 
     - You can feed it 2D or 3D data. It will always average over angles to make 
    an isotropic spectrum
     - This version does the 'classic' calculation of the 1d spectrum, not taking into account the fact that for large |k| less modes are represented.
    """
    # Construct the modes
    ns = []
    for i in range(domain.dim):
        k = domain.elements(i)
        L = np.sum(domain.bases[i].grid_spacing())
        n = np.round(k*(L/(2*np.pi))).astype(np.int)
        ns.append(n)
        
    n2 = np.sqrt(np.sum([n**2 for n in ns],axis=0))
     
    # Find maximum mode number modulus:
    n2max = np.round(np.amax(n2)).astype(np.int)
    n2max_all = 0.0
    n2max_all = comm.allreduce(n2max,op=MPI.MAX)
              
    # square them and add them up to get energy density        
    KEs = 0.5*np.sum([np.abs(field['c'])**2 for field in fields],axis=0)
    KEs *= (2 - (ns[0] == 0))  # Double-count positive kx because of R2C transform

    # Construct arrays for final output
    modn = np.arange(1,n2max_all+1)
    En = np.zeros(n2max_all)   
    
    for k in modn:
        En[modn==k] = np.sum(KEs[(np.round(n2).astype(np.int))==k])
        
    # Do zeroth first:
    En[modn==1] += np.sum(KEs[(np.round(n2).astype(np.int))==0])
            
    Entot = comm.allreduce(En,op=MPI.SUM)

    return [modn,Entot]

def spec1d(domain,fields,modes=True):
    """
    [modn,En] = spec1d(domain,[u,v,w]) or
    [modn,En] = spec1d(domain,[u,v])

    input: domain, fields (list of dedalus fields, like [u,v,w])
    output: [modn,En], 1-dimensional power spectrum of u['c']**2+v['c']**2 + ...
    
    Options:
    modes (=True by default): if True, constructs the 1d spectrum based on mode numbers.
    otherwise, it constructs the 1d spectrum based on wave numbers.
    
    NOTES/WARNINGS: 
     - You can feed it 2D or 3D data. It will always average over angles to make 
    an isotropic spectrum
     - If 'modes=False', i.e. you're using wavenumbers, it won't work because this assumes isotropic domain! If you want this to work, you will have to modify the code slightly :)
     - This version compensates the fact that some modes are less represented in a spherical shell with a cartesian lattice of wavenumbers. 
     - This version is also 10x faster than the spec1d_classic, so use this one if possible!
    """
    if modes:
        # Construct the modes
        ns = []
        for i in range(domain.dim):
            k = domain.elements(i)
            L = np.sum(domain.bases[i].grid_spacing())
            n = np.round(k*(L/(2*np.pi))).astype(np.int)
            ns.append(n)
            dn = 1
    else:
        # Construct the wavenumbers
        ns = []
        for i in range(domain.dim):
            k = domain.elements(i)
            L = np.sum(domain.bases[i].grid_spacing())
            ns.append(k)
            dn = 2 * np.pi / L
        
    n2 = np.sqrt(np.sum([n**2 for n in ns],axis=0)) # not n^2, just |n|
    
    # Find maximum mode number modulus:
    n2max = np.round(np.amax(n2)).astype(np.int)
    n2max_all = comm.allreduce(n2max,op=MPI.MAX)
    
    # Take histogram of modes lying in each spherical shell
    # Make the first interval from 0 to 1.5, then the rest add one to it.
    bins = np.concatenate(([0.0],np.arange(1.5, n2max_all+1.5, 1)))
        
    hist_samples, _ = np.histogram(n2, bins=bins)
    hist_samples = comm.allreduce(hist_samples,op=MPI.SUM)
              
    # square feilds and add them up to get energy density        
    KEs = 0.5*np.sum([np.abs(field['c'])**2 for field in fields],axis=0) / dn**(domain.dim)
    if domain.dim==2:
        KEs_1d = KEs * np.pi * n2
    elif domain.dim==3:
        KEs_1d = KEs * 4 * np.pi * n2**2
    else:
        print("ERROR: can only compute power spectrum from 2 or 3 dimensional data.")
    KEs_1d *= (2 - (ns[0] == 0))  # Double-count positive kx because of R2C transform

    # Plot histogram
    pow_samples, _ = np.histogram(n2, bins=bins, weights=KEs_1d)
    pow_samples = comm.allreduce(pow_samples,op=MPI.SUM)
    
    Entot = pow_samples / hist_samples
            
    return [np.arange(1,len(Entot)+1),Entot]
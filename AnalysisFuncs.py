import numpy as np
import matplotlib.pyplot as plt

def CorrThetaF(n, ThetaF):
    # biais correction for thetaF see Xu and Fu 2004
    if ThetaF>15:
        A = (1.1675 + 3.3232/n + 63.698/(n**2))
        B = 0.2569
    else:
        A = (1.1313 + 3.4882/n + 28.2878/(n**2))
        B = 0.3998
    UnbiasedThetaF = ((np.sqrt(B**2 + 4*A*ThetaF) - B)/(2*A))**2
    return UnbiasedThetaF

def Unique_genotype(N, nuc, loci_to_compare = None):
    """return the proportion and the eveness of unique genotypes calculated using the loci index given in loci to compare if none are provided uses all microsat loci [2:]"""
    # Array with 1 genotype = 1 line
    if loci_to_compare==None:
        loci_to_compare = [i for i in range(2, len(nuc[0,:]))]
    nloci = len(loci_to_compare)
    gen_array = nuc[:,loci_to_compare].reshape(N, 2*nloci)
    uniques, counts = np.unique(gen_array, axis=0, return_counts=True)

    # Number of genotypes
    G = len(uniques)
    # unique proportion
    P = np.sum(counts==1)/N
    # eveness with the simpson 1 / D
    E_1D = (1/np.sum((counts/N)**2)) / G
    return G, P, E_1D

def NAll(nuc, mean_per_locus = True):
    """return number of allele for each locus"""
    nloc = len(nuc[0,])
    NbAll = [None]*(nloc-2)
    for idxlocus in range(2, nloc):
        NbAll[idxlocus-2] = len(np.unique(nuc[:,idxlocus]))
    if mean_per_locus:
        NbAll = np.mean(NbAll)
    return NbAll

def Rst_Sexes(N, nuc, mean_per_locus = True):
    Nm = np.sum(nuc[:,0])
    Nf = N-Nm
    # subpop of females and males
    Idx_m = (nuc[:, 0]==1).nonzero()[0]
    Idx_m = np.concat([Idx_m, Idx_m-1])
    Idx_f = np.ones(N*2, dtype=bool)
    Idx_f[Idx_m] = False

    S = np.var(nuc[:, 2:], axis=0)
    Sm = np.var(nuc[Idx_m, 2:], axis=0)
    Sf = np.var(nuc[Idx_f, 2:], axis=0)
    Sw = Sm+Sf/2
    Rst = 1 - Sw/S
    if mean_per_locus:
        Rst = np.mean(Rst)
    return Rst

def heterozigoty(N, nuc, mean_per_locus = True):
    """return the expected heterozigoty at HW equilibrium"""
    nloc = len(nuc[0,])
    H = [None]*(nloc-2)
    for idxlocus in range(2, nloc):
        freqs = np.unique(nuc[:,idxlocus], return_counts=True)[1]/(N*2)
        H[idxlocus-2] = 1-np.sum(freqs**2)
    if mean_per_locus:
        H = np.mean(H)
    return H

def Fst_Sexes(N, nuc, mean_per_locus = True):
    Nm = np.sum(nuc[:,0])
    Nf = N-Nm
    # subpop of females and males
    Idx_m = (nuc[:, 0]==1).nonzero()[0]
    Idx_m = np.concat([Idx_m, Idx_m-1])
    Idx_f = np.ones(N*2, dtype=bool)
    Idx_f[Idx_m] = False

    Hm = np.array(heterozigoty(Nm, nuc[Idx_m], False))
    Hf = np.array(heterozigoty(Nf, nuc[Idx_f], False))
    Ht = np.array(heterozigoty(N, nuc, False))
    Hs = (Nm*Hm+Nf*Hf)/N
    Fst = 1 - Hs/Ht
    if mean_per_locus:
        Fst = np.mean(Fst)
    return Fst

def SexNumber(Model, N, nuc, cyt = None):
    """return [Nm, Nf] for dioecious models and [Nh, Nm, Nmcms, Nf] for trioecy"""
    if Model=="M1" or Model=="M2":
        Nm = np.sum(nuc[:,0])
        Nf = N-Nm
        return [Nm, Nf]
    elif Model=="M3":
        Nm,Nf,Nmcms,Nh = 0,0,0,0
        for i in range(N):
            if nuc[i*2+1,0]==1:
                if cyt[i,0]==1:
                    Nmcms+=1
                Nm+=1
            if nuc[i*2+1,0]==0 and cyt[i,0]==1:
                Nf+=1
            if nuc[i*2+1,0]==0 and cyt[i,0]==0:
                Nh+=1
        return  [Nh, Nm, Nmcms, Nf]

def NucEffectiveSize(N, Model, nuc, Musat, Muiloc, nMsats, cyt = None):
    NeEstimations = [] # [SR, Inf Locus, thetaF, thetaV]
    if Model == "M1" or Model == "M2":
        Nm, Nf = SexNumber(Model, N, nuc)
        NeEstimations.append(4*Nm*Nf/(N))
    elif Model == "M3":
        NeEstimations.append(None)
    freqs = np.unique(nuc[:, 1], return_counts=True)[1]/(N*2)
    He = 1 - np.sum(freqs**2)
    theta = He/(1-He)
    NeEstimations.append(theta/(4*Muiloc))
    NeF = [None] * nMsats
    NeV = [None] * nMsats
    for i in range(nMsats):
        freqs = np.unique(nuc[:, i+2], return_counts=True)[1]/(N*2)
        F = np.sum(freqs**2)
        thetaF = 0.5*(1/(F**2) - 1)
        NeF[i] = CorrThetaF(N, thetaF)/(4*Musat)
        thetaV = np.var(nuc[:, i+2])*2
        NeV[i] = thetaV/(4*Musat)
    NeEstimations.append(np.mean(NeF))
    NeEstimations.append(np.mean(NeV))
    return NeEstimations



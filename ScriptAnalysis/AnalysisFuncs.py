import numpy as np
import matplotlib.pyplot as plt

def load_sim(dirname, N):
    nuc = np.load(dirname+"nuc"+str(6*N)+".npy")
    cyt = np.load(dirname+"cyt"+str(6*N)+".npy")
    return nuc, cyt

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
    if np.any(loci_to_compare==None):
        loci_to_compare = [i for i in range(2, len(nuc[0,:]))]
    nloci = len(loci_to_compare)

    #Transform nuc into ordoned ensemble of alleles pair
    gen_array = np.zeros((N, 2*nloci), dtype=int)
    for i in range(nloci):
        gen_array[:,i*2:i*2+2] = nuc[:, 2+i].reshape((N, 2))

    #ordination
    gen_array = np.sort(gen_array.reshape(N, nloci, 2), axis=2).reshape(N, 2*nloci)
    uniques, counts = np.unique(gen_array, axis=0, return_counts=True)

    # Number of genotypes
    G = len(counts)
    # unique proportion
    Pu = np.sum(counts==1)/G
    Pi = G/N
    # eveness with the simpson 1 / D
    E_1D = (1/np.sum((counts/N)**2)) / G
    return Pu, Pi, E_1D

def split_nuc(nuc, trioecy = False, cyt = None):
    males_nuc = np.repeat(nuc[1::2,0]==1, 2)
    if trioecy:
        cms = np.repeat(cyt[:,0]==1, 2) #nuc position of individual with cms
        males_cms = np.logical_and(males_nuc, cms)
        males_Ncms = np.logical_and(males_nuc, np.bool_(1-cms))
        females = np.logical_and(np.bool_(1-males_nuc), cms)
        hermaphrodites = np.logical_and(np.bool_(1-males_nuc), np.bool_(1-cms))
        return nuc[hermaphrodites,:], nuc[males_Ncms, :],  nuc[males_cms, :],  nuc[females, :]
    else:
        return nuc[males_nuc,:], nuc[np.bool_(1-males_nuc), :]

def NAll(nuc, mean_per_locus = True, per_sex_class = [False], Ncor = None, repeats = 100):
    """return number of allele for each locus
    per_sex_class = [True/False -> calculate per class, True/False -> Trioecy, cyt]"""
    nloc = nuc.shape[1]
    if Ncor!=None and Ncor!= "no" and Ncor<=5:
        return np.nan
    NbAll = np.zeros(nloc-2, dtype=float)
    if per_sex_class[0]:
        if per_sex_class[1]:
            s_nuc = split_nuc(nuc, True, per_sex_class[2])
        else:
            s_nuc = split_nuc(nuc)
        min_size = s_nuc[0].shape[0]
        for split in s_nuc[1:]:
            if min_size> split.shape[0]:
                min_size = split.shape[0]
        list_Nall = []
        for split in s_nuc:
            Cor_size = min_size//2
            if Ncor=="no": Cor_size = None
            list_Nall.append(NAll(split, mean_per_locus, Ncor = Cor_size, repeats = repeats))
        return list_Nall
    else:

        for idxlocus in range(2, nloc):
            if Ncor==None:
                NbAll[idxlocus-2] = len(np.unique(nuc[:,idxlocus]))
            #correction for sample size
            else:
                sum_loc = 0
                N_ind = nuc.shape[0]//2
                for i in range(repeats):
                    idxindividuals = np.zeros(Ncor*2, int)
                    idxindividuals[0::2] = np.random.choice(N_ind, replace=True, size=Ncor)*2
                    idxindividuals[1::2] = idxindividuals[0::2]+1
                    sum_loc += len(np.unique(nuc[idxindividuals,idxlocus]))
                NbAll[idxlocus-2] = sum_loc/repeats

        if mean_per_locus:
            NbAll = np.mean(NbAll)
    return NbAll

def heterozigoty(N, nuc, mean_per_locus = True, loci_iterator=None):
    """return the expected heterozigoty at HW equilibrium"""
    nloc = len(nuc[0,])
    H = np.zeros(nloc-2, dtype=float)
    if loci_iterator==None: loci_iterator = range(2, nloc)
    for idxlocus in loci_iterator:
        freqs = np.unique(nuc[:,idxlocus], return_counts=True)[1]/(N*2)
        H[idxlocus-2] = 1-np.sum(freqs**2)
    if mean_per_locus:
        H = np.mean(H)
    return H

def Fis(N, nuc, global_fis = True):
    nloc = len(nuc[0,])
    Fis = np.zeros((nloc-2), dtype=float)# only microsats
    for ind in range(N):
        Fis += np.int8(nuc[ind*2,2:] != nuc[ind*2+1,2:]) # count the number of heterozygous individuals
    Fis = Fis/N
    He = heterozigoty(N, nuc, mean_per_locus = False)
    if global_fis:
        Fis = Fis[He!=0]
        Fis = 1 - np.sum(Fis)/np.sum(He[He!=0])
    else:
        Fis[He!=0] =  1 - Fis[He!=0]/He[He!=0]
        Fis[He==0] = np.nan
    return Fis


def SexNumber(Model, N, nuc, cyt = None):
    """return [Nm, Nf] for dioecious models and [Nh, Nm, Nmcms, Nf] for trioecy"""
    if Model=="M1" or Model=="M2":
        Nm = np.sum(nuc[:,0])
        Nf = N-Nm
        return [Nm, Nf]
    elif Model=="M3":
        H, M, C, F = split_nuc(nuc, True, cyt)
        Nh = H.shape[0]
        Nm = M.shape[0]
        Nmcms = C.shape[0]
        Nf = F.shape[0]
        return  np.array([Nh, Nm, Nmcms, Nf])//2

def NucEffectiveSize(N, Model, nuc, Musat, Muiloc, nMsats, cyt = None):
    NeEstimations = [] # [thetaF, thetaV]
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

def indID_to_nucID(indID):
    nuc_ID = np.zeros(indID.shape[0]*2, int)
    nuc_ID[::2] = indID*2
    nuc_ID[1::2] = nuc_ID[::2] + 1
    return nuc_ID

def HtHs(N_i, N_j, nuc_i, nuc_j, nMsats=10):
    """return Ht and Hs see Nei 1973"""
    w_i, w_j = N_i/(N_i+N_j), N_j/(N_i+N_j) # sexe classes weights
    Jt, Js = np.zeros(nMsats, float), np.zeros(nMsats, float) #Homozygoty / probability identity
    pooled_nuc = np.concat((nuc_i, nuc_j))
    for sat_ID in range(2, nMsats+2):
        Alleles_array = np.unique(pooled_nuc[:, sat_ID])
        x_i = np.sum(nuc_i[:, sat_ID][:,None] == Alleles_array, axis=0)/(N_i*2)
        x_j = np.sum(nuc_j[:, sat_ID][:,None] == Alleles_array, axis=0)/(N_j*2)
        Jt[sat_ID-2] = np.sum((w_i*x_i + w_j*x_j)**2)
        Js[sat_ID-2] = w_i*np.sum(x_i**2) + w_j*np.sum(x_j**2)
    Ht, Hs = np.mean(1-Jt), np.mean(1-Js) # Sum on all loci cf Nei 1977
    return Ht, Hs

def Dist_mf(nuc, trioecy, cyt=None):
    snuc = split_nuc(nuc, trioecy, cyt)
    N  = [snuc[i].shape[0]//2 for i in range(len(snuc))]
    if trioecy:
        morph_m, morph_f = np.concat(snuc[1:3]), snuc[3]
        Nm, Nf = np.sum(N[0:3]), N[3]
        Ht, Hs = HtHs(Nm, Nf, morph_m, morph_f)
    else:
        Ht, Hs = HtHs(N[0], N[1], snuc[0], snuc[1])
    return [1 - Hs/Ht, 2 * (Ht-Hs)/(1 - Hs)]

def Dist_sexmatrix(nuc, cyt, test = False, n = 100, nMsats = 10):
    """return the Dist matrix using lower triangle Nei 1973 (CalcGst) and upper triangle Jost 2008 (CalcJostD) methods
      H M C F
    H
    M
    C
    F
    """
    mat = np.zeros((4, 4))
    snuc = split_nuc(nuc, True, cyt) #[h, m, c, f]
    N = snuc[0].shape[0]//2, snuc[1].shape[0]//2, snuc[2].shape[0]//2, snuc[3].shape[0]//2 #[Nh, Nm, Nc, Nf]
    for i in range(3):
        for j in range(i+1, 4):
            Ht,Hs = HtHs(N[i], N[j], snuc[i], snuc[j], nMsats)
            mat[j, i] = 1 - Hs/Ht #Gst
            mat[i, j] = 2 * (Ht-Hs)/(1 - Hs) #Jost's D

    if test:
        Ntot = np.sum(N)
        distrmat = np.zeros((n,4,4), float)
        for i in range(3):
            for j in range(i+1, 4):
                pooled_nuc = np.concat((snuc[i], snuc[j]))
                for rep in range(n):
                    perm = np.random.permutation(N[i] + N[j])
                    IDi = indID_to_nucID(perm[:N[i]])
                    IDj = indID_to_nucID(perm[N[i]:])
                    Ht,Hs = HtHs(N[i], N[j], pooled_nuc[IDi,], pooled_nuc[IDj,], nMsats)
                    distrmat[rep, j, i] = 1 - Hs/Ht #Gst
                    distrmat[rep, i, j] = 2 * (Ht-Hs)/(1 - Hs) #Jost's D
        return mat, distrmat
    return mat

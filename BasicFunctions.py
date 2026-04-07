import numpy as np


def mutations_iloc(N, Mu, pop, mut_ID, NucCyt):
    nMu = np.random.binomial(N, Mu)
    idx = np.random.choice(N, nMu, replace = False)
    old = mut_ID[NucCyt]
    mut_ID[NucCyt] += nMu
    pop[idx] = np.arange(old, mut_ID[NucCyt])
    return

def mutation_msat(N2, nMsats, Mu, msat_array):
    TotMsats = (N2*nMsats)
    nMut = np.random.binomial(TotMsats, Mu)
    MutIndex = np.random.choice(TotMsats, nMut, replace = False)
    MutIndex = np.unravel_index(MutIndex, (N2, nMsats))
    Additional_repeats = np.random.choice([-1, 1], nMut, replace = True)
    msat_array[MutIndex] += Additional_repeats
    return

def recombi_ind(c1c2, nMsats):
    """return recombination of c1 with c2 whith nCo crossing over
    c1c2: array 2*(nMsats+2) with the two chromosomes"""
    idx = np.random.randint(2, size = nMsats+2)
    recombined = np.copy(c1c2[idx, range(nMsats+2)])
    return recombined

def Clonality(Nsurv, c, tmp_nuc, tmp_cyt):
    Ncopy = np.random.binomial(Nsurv,c)
    IdxGenets = np.random.choice(Nsurv, Ncopy, replace = False)
    Ncopy2 = Ncopy*2
    IdxgenetsNuc = np.zeros(Ncopy2, dtype=int)
    IdxgenetsNuc[0:Ncopy2:2] = IdxGenets*2
    IdxgenetsNuc[1:Ncopy2:2] = IdxgenetsNuc[0:Ncopy2:2]+1
    tmp_nuc[Nsurv2:(Nsurv2+Ncopy2), : ] = np.copy(tmp_nuc[IdxgenetsNuc, : ])
    tmp_cyt[Nsurv:(Nsurv+Ncopy)] = np.copy(tmp_cyt[IdxGenets])
    return Ncopy

def ClonalityPerenity(param, N, nuc, cyt, tmp_nuc, tmp_cyt, Age):
    """First draw the survivants between two years (no biais) and then produce clonal copy.
return the index after the last individual added"""
    # Drawing surviving individuals
    MaxAge, p, c = param
    #age Weights
    Weights = np.zeros(N) + 1
    Weights[Age==MaxAge] = 0
    NCandidates = N - np.sum(Age==MaxAge)

    Nsurv = np.random.binomial(NCandidates*p)
    IdxsurvInd = np.random.choice(N, Nsurv, replace = False,  p = Weights/np.sum(Weights))
    Nsurv2 = Nsurv*2
    IdxsurvNuc = np.zeros(Nsurv2, dtype=int)
    IdxsurvNuc[0:Nsurv2:2] = IdxsurvInd*2
    IdxsurvNuc[1:Nsurv2:2] = IdxsurvNuc[0:Nsurv2:2]+1
    tmp_nuc[0:Nsurv2, : ] = np.copy(nuc[IdxsurvNuc, : ])
    tmp_cyt[0:Nsurv] = np.copy(cyt[0:Nsurv])
    tmp_Age = np.zeros(N, dtype=int)
    tmp_Age[IdxsurvInd] = Age[IdxsurvInd]+1
    Age[:] = np.copy(tmp_Age)

    # Drawing the survivants producing separate ramets
    Ncopy = Clonality(Nsurv, c, tmp_nuc, tmp_cyt)
    return Nsurv+Ncopy

def ClonalityPerenityBiased(param, N, nuc, cyt, tmp_nuc, tmp_cyt, Age):
    """First draw the survivants between two years (no biais) and then produce clonal copy.
return the index after the last individual added"""
    # Drawing surviving individuals
    MaxAge, p, c, pm = param

    #Weights calcul
    Weights = np.zeros(N) + 1
    old = [Age==MaxAge]
    Weights[old] = 0

    males = np.bool(nuc[1::2,0]) # Array True for male & False for female
    females = np.bool(1 - maleIDs)
    # candidates for each sex
    Nm = np.sum(maleIDs)
    Nf = N - Nm
    Nmc = Nm - np.sum(old[maleIDs])
    Nfc = Nf - np.sum(old[femaleIDs])

    NsurvM = np.random.binomial(Nmc*pm)
    NsurvF = np.random.binomial(Nfc*p)
    IdxsurvInd = np.concat(np.random.choice(np.arange(N)[males], NsurvM, replace = False,  p = Weights[males]/np.sum(Weights[males])),
                           np.random.choice(np.arange(N)[females], NsurvF, replace = False,  p = Weights[females]/np.sum(Weights[females])))

    # Drawing the survivants producing separate ramets
    Nsurv = NsurvF + NsurvM
    Ncopy = Clonality(Nsurv, c, tmp_nuc, tmp_cyt)
    return Nsurv+Ncopy

def No_CP(param, N, nuc, cyt, tmp_nuc, tmp_cyt, Age):
    return 0

# def recombi(Pmsats, c1, c2, nCo, nMsats): #pmsats -> microsat position
#     """return recombination of c1 with c2 whith nCo crossing over"""
#     recombined = np.zeros(nMsats+1, dtype=int)
#     recombined[1:] = c1[:] # np.copy
#     CoPos = np.sort(np.random.rand(nCo)) #crossing over positions
#     j = 0
#     for i in range(nCo):
#         while(j<nMsats and Pmsats[j]<=CoPos[i]):
#             if i%2==1:
#                 recombined[j+1] = c2[j]
#             j+=1
#     if nCo%2==1: # TESTER
#         for i in range(j, nMsats):
#             recombined[i+1] = c2[i]
#     return recombined

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

def ClonalityUB(Nsurv, Nsurv2, param, tmp_nuc, tmp_cyt):
    c = param[2]
    Ncopy = np.random.binomial(Nsurv, c)
    IdxGenets = np.random.choice(Nsurv, Ncopy, replace = False)
    Ncopy2 = Ncopy*2
    IdxgenetsNuc = np.zeros(Ncopy2, dtype=int)
    IdxgenetsNuc[0:Ncopy2:2] = IdxGenets*2
    IdxgenetsNuc[1:Ncopy2:2] = IdxgenetsNuc[0:Ncopy2:2]+1
    tmp_nuc[Nsurv2:(Nsurv2+Ncopy2), : ] = np.copy(tmp_nuc[IdxgenetsNuc, : ])
    tmp_cyt[Nsurv:(Nsurv+Ncopy)] = np.copy(tmp_cyt[IdxGenets])
    return Ncopy

def ClonalityB(Nsurv, Nsurv2, param, tmp_nuc, tmp_cyt):
    c, cm = param[2:]
    males = np.bool(tmp_nuc[1:Nsurv2:2, 0])
    females = np.bool(1-males)

    Nsurvm = np.sum(males)
    Nsurvf = Nsurv-Nsurvm
    Ncopym = np.random.binomial(Nsurvm, cm)
    Ncopyf = np.random.binomial(Nsurvf, c)

    IdxGenets = np.concat((np.random.choice(np.arange(Nsurv)[males], Ncopym, replace = False),
                           (np.random.choice(np.arange(Nsurv)[females], Ncopyf, replace = False))))
    Ncopy = Ncopym+Ncopyf
    Ncopy2 = Ncopy*2
    IdxgenetsNuc = np.zeros(Ncopy2, dtype=int)
    IdxgenetsNuc[0:Ncopy2:2] = IdxGenets*2
    IdxgenetsNuc[1:Ncopy2:2] = IdxgenetsNuc[0:Ncopy2:2]+1
    tmp_nuc[Nsurv2:(Nsurv2+Ncopy2), : ] = np.copy(tmp_nuc[IdxgenetsNuc, : ])
    tmp_cyt[Nsurv:(Nsurv+Ncopy)] = np.copy(tmp_cyt[IdxGenets])
    return Ncopy

def ClonalityPerenity(param, N, nuc, cyt, tmp_nuc, tmp_cyt, Age, Clfunc):
    """First draw the survivants between two years (no biais) and then produce clonal copy.
return the index after the last individual added"""
    # Drawing surviving individuals
    MaxAge, p = param[:2]
    #age Weights
    Weights = np.zeros(N) + 1
    Weights[Age==MaxAge] = 0
    NCandidates = N - np.sum(Age==MaxAge)

    Nsurv = np.random.binomial(NCandidates, p)
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
    Ncopy = Clfunc(Nsurv, Nsurv2, param, tmp_nuc, tmp_cyt)
    return Nsurv+Ncopy

def No_CP(param, N, nuc, cyt, tmp_nuc, tmp_cyt, Age, Clfunc):
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

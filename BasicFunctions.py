import numpy as np

def mutations_iloc(N, pop, Mu):
    nMu = np.random.poisson(N*Mu)
    idx = np.random.choice(range(N), nMu, replace = True)
    pop[idx] = np.random.rand(nMu)*1e16 #use int for practicity
    return pop

def mutation_msat(seq, nMut, nMsats):
    MutMsatsIndex = np.random.choice(range(nMsats), nMut, replace = True) + 2 #+2 for the sex locus and the neutral locus
    Additional_repeats = np.random.choice([-1, 1], nMut, replace = True)
    seq[MutMsatsIndex] = seq[MutMsatsIndex] + Additional_repeats
    return

def recombi_ind(c1c2, nMsats):
    """return recombination of c1 with c2 whith nCo crossing over
    c1c2: array 2*(nMsats+2) with the two chromosomes"""
    #optimisable
    recombined = np.zeros(nMsats+2, dtype=int)
    idx = np.random.randint(2, size = nMsats+2)
    recombined = c1c2[idx, range(nMsats+2)]
    return recombined

# def recombi(Pmsats, c1, c2, nCo, nMsats): #pmsats -> microsat position
#     """return recombination of c1 with c2 whith nCo crossing over"""
#     recombined = np.zeros(nMsats+1, dtype=int)
#     recombined[1:] = c1[:]
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

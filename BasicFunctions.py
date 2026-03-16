import numpy as np

def mutations_cyt(N, pop_cyt, Mucyt):
    nMu = np.random.poisson(N*Mucyt)
    idx = np.random.choice(range(N), nMu)
    pop_cyt[idx] = np.random.rand(nMu)
    return

def mutation_msat(seq, nMut, nMsats):
    MutMsatsIndex = np.random.choice(range(nMsats), nMut) + 1 #+1 for the sex locus
    Additional_repeats = np.random.choice([-1, 1], nMut)
    seq[MutMsatsIndex] = seq[MutMsatsIndex] + Additional_repeats
    return

def recombi_ind(c1c2, nMsats):
    """return recombination of c1 with c2 whith nCo crossing over
    c1c2: array 2*(nMsats+1) with the two chromosomes"""
    #optimisable
    recombined = np.zeros(nMsats+1, dtype=int)
    idx = np.random.randint(2, size = nMsats+1)
    recombined = c1c2[idx, range(nMsats+1)]
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

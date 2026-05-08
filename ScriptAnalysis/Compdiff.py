import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *

from AnalysisFuncs import *

N=1000
#trio
#1
D_fm_trio1p = []
for i in range(100):
    dirname = "[...]"+"/M3_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)
    H = SexNumber("M3", N, nuc, cyt = cyt)[0]
    if H!=0:
        D_fm_trio1p.append(Dist_mf(nuc, True, cyt)[1])

#2
D_fm_trio2p = []
for i in range(100):
    dirname ="[...]"+"/M3_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)
    H = SexNumber("M3", N, nuc, cyt = cyt)[0]
    if H!=0:
        D_fm_trio2p.append(Dist_mf(nuc, True, cyt)[1])

# Null
D_fm_null = []
for i in range(100):
    dirname ="[...]"+"/M1_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)
    D_fm_null.append(Dist_mf(nuc, False)[1])

# Cpb
D_fm_cpb = []
for i in range(100):
    dirname ="[...]"+"/M1_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)
    D_fm_cpb.append(Dist_mf(nuc, False)[1])

#Sm 0.5
D_fm_sm05 = []
for i in range(100):
    dirname ="[...]"+"/M2_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)
    D_fm_sm05.append(Dist_mf(nuc, False)[1])

# Sm 0.75
D_fm_sm075 = []
for i in range(100):
    dirname ="[...]"+"/M2_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)
    D_fm_sm075.append(Dist_mf(nuc, False)[1])


def affichage(arr,comp_to = D_fm_null):
    w = mannwhitneyu(arr, comp_to)
    print("(w =", w.statistic, "p =", str(w.pvalue)+")")

affichage(D_fm_trio1p)
affichage(D_fm_trio2p)
affichage(D_fm_cpb)
affichage(D_fm_sm05)
affichage(D_fm_sm075)

plt.hist(D_fm_cpb)
plt.hist(D_fm_null)

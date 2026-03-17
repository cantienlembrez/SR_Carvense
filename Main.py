import numpy as np
from os import mkdir
from datetime import datetime

from BasicFunctions import *
from Models import *

def CorrThetaF(n, ThetaF):
    # biais correction see (Xu and Fu 2004)
    if ThetaF>15:
        A = (1.1675 + 3.3232/n + 63.698/(n**2))
        B = 0.2569
    else:
        A = (1.1313 + 3.4882/n + 28.2878/(n**2))
        B = 0.3998
    UnbiasedThetaF = ((np.sqrt(B**2 + 4*A*ThetaF) - B)/(2*A))**2
    return UnbiasedThetaF

def run(N, IDmsats, Musat, Muiloc, Gmax, Sm=0.5):
    """N: number of individuals,
    IDmsats: scaled position of microsatellites loci on chromosome (increasing order)
    Musat: microsatellites mutation rate
    Muiloc: mutation rate for infinite allele loci
    Gmax: generation ending simulation
    Sm: survival rate of male in model 2
    """
    N2 = N*2
    nMsats = len(IDmsats)

    #initializing microsatellites each individual share the same haplotype at start
    #each individual is two line one for each chromosome first mother herited haplotype
    #the first column contain the information for the sex determination locus 0 : female recessive
    # and the second contain an neutral locus (infinite allele locus)
    nuc = np.zeros((N2, nMsats+2),dtype=int)
    #infinite allele model one neutral locus
    cyt = np.zeros((N),dtype=int)

    tmp_nuc = np.zeros((N2, nMsats+2),dtype=int)
    tmp_cyt = np.zeros((N),dtype=int)
    repeats = np.random.poisson(20, nMsats) + 5
    for chr_id in range(N2):
        nuc[chr_id,2:] = repeats
        #initial sex ratio of 0.5
        if chr_id<N and chr_id%2==1:
            nuc[chr_id,0] = 1

    f = open("Simulations/"+datetime.today().strftime('%Y-%m-%d_%H:%M') +".txt", "a")
    #writing parameters and initial state
    params =[N, IDmsats, Musat, Muiloc, Gmax, Sm]
    for i in  range(len(params)):
        f.write(run.__code__.co_varnames[i]+":"+str(params[i])+"    ")
    f.write("\n")
    f.write("      " + "   ".join(["_".join([str(i), str(IDmsats[i])],)  for i in range(nMsats)])+"\n")
    f.write("Anc   "+"   ".join([str(r) for r in repeats])+"\n")

    for g in range(1, Gmax+1):
        #gen_incr(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt)
        gen_incr_m2(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt, Sm)
        nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)

    Nm = np.sum(nuc[1:N2:2, 0])
    Nf = N-Nm
    #effective size estimation for nuclear loci
    lthetaF =[None]*nMsats
    lthetaV =[None]*nMsats
    #inf allele
    freqs = np.unique(nuc[:, 1], return_counts=True)[1]/N2
    He = 1 - sum(freqs**2)
    theta = He/((1-He))
    print(He)
    print("theta from infinite allele locus :" +str(theta))
    print("estimated effective pop size : " + str(theta/(4*Muiloc)))
    for i in range(nMsats):
        f.write(IDmsats[i]+":\n")
        tab_sat = np.unique(nuc[:, i+2], return_counts=True)
        f.write("allele"+"  "+"counts\n")
        for j in range(len(tab_sat[0])):
            f.write(str(tab_sat[0][j])+"    "+str(tab_sat[1][j])+"\n")
        F = sum((tab_sat[1]/N2)**2)
        thetaF = 0.5*(1/(F**2) - 1)
        lthetaF[i] = thetaF
        Vs = np.var(nuc[:, i+2])
        thetaV = 2*Vs
        lthetaV[i] = thetaV

    thetaF = np.mean(lthetaF)
    thetaV = np.mean(lthetaV)
    print("mean thetaF :", thetaF)
    print("corrected thetaF :", CorrThetaF(N2,thetaF))
    print("estimated pop size", CorrThetaF(N2, thetaF)/(4*Musat))
    print("mean thetaV :", thetaV)
    print("estimated pop size", thetaV/(4*Musat))
    print("sex ratio :" +str(Nm/(Nm+Nf)))
    Ne = (Nm*Nf*4)/(Nm+Nf)
    print("population size estimated from sex ratio : " + str(Ne))

    #effective size estimation for cytoplasmic loci
    alleles, freqs  = np.unique(cyt, return_counts=True)
    freqs=freqs/N
    He = 1 - sum(freqs**2)
    thetaC = He/((1-He))
    print("cyt theta :" +str(thetaC))
    print("estimated effective pop size : " + str(thetaC/(2*Muiloc)))
    print("number of female:",str(Nf))
    print("estimated pop size:",str(Nf/2))
    f.close()
    return

try:
    mkdir("Simulations")
except:
    pass

Musat = 1e-2
N = 2000
print("theta theo : " + str(4*N*Musat))
run(N,[chr(i+65) for i in range(25)], Musat, Musat, 9000)

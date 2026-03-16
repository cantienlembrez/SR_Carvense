#!/bin/bash

import numpy as np
from os import mkdir
from datetime import datetime

from BasicFunctions import *
from Models import *



def run(N, IDmsats, Musat, Mucyt, Gmax, Sm=0.5):
    """N: number of individuals,
    IDmsats: scaled position of microsatellites loci on chromosome (increasing order)
    Musat: microsatellites mutation rate
    Mucyt: cytoplasm mutation rate
    Gmax: generation ending simulation
    Sm: survival rate of male in model 2
    """
    N2 = N*2
    nMsats = len(IDmsats)

    #initializing microsatellites each individual share the same haplotype at start
    #each individual is two line one for each chromosome first mother herited haplotype
    #the first column contain the information for the sex determination locus 0 : female recessive
    # and the second contain an neutral locus (infinite allele locus)
    nuc = np.zeros((N2, nMsats+1),dtype=int)
    #infinite allele model one neutral locus
    cyt = np.zeros((N),dtype=float)

    tmp_nuc = np.zeros((N2, nMsats+1),dtype=int)
    tmp_cyt = np.zeros((N),dtype=float)
    repeats = np.random.poisson(20, nMsats) + 5
    for chr_id in range(N2):
        nuc[chr_id,1:] = repeats

        #initial sex ratio of 0.5
        if chr_id<N and chr_id%2==1:
            nuc[chr_id,0] = 1

    f = open("Simulations/"+datetime.today().strftime('%Y-%m-%d_%H:%M') +".txt", "a")
    #writing parameters and initial state
    params =[N, IDmsats, Musat, Mucyt, Gmax, Sm]
    for i in  range(len(params)):
        f.write(run.__code__.co_varnames[i]+":"+str(params[i])+"    ")
    f.write("\n")
    f.write("      " + "   ".join(["_".join([str(i), str(IDmsats[i])],)  for i in range(nMsats)])+"\n")
    f.write("Anc   "+"   ".join([str(r) for r in repeats])+"\n")


    for g in range(1, Gmax+1):
        gen_incr(N, N2, nMsats, Musat, Mucyt, nuc, cyt, tmp_nuc, tmp_cyt)
        #gen_incr_m2(N, N2, IDmsats, nMsats, Musat, Mucyt, nuc, cyt, tmp_nuc, tmp_cyt, Sm)
        nuc, cyt = tmp_nuc, tmp_cyt
        if g%1000==0:
            print(g)


    Nm = np.sum(nuc[1:N2:2, 0])
    Nf = N-Nm
    #effective size estimation for nuclear loci
    lthetaF =[None]*nMsats
    lthetaV =[None]*nMsats
    for i in range(nMsats):
        f.write(IDmsats[i]+":\n")
        tab_sat = np.unique(nuc[:, i+1], return_counts=True)
        f.write("allele"+"  "+"counts\n")
        for j in range(len(tab_sat[0])):
            f.write(str(tab_sat[0][j])+"    "+str(tab_sat[1][j])+"\n")
        print("microsatellite locus "+ IDmsats[i] + ": ")
        F = sum((tab_sat[1]/N2)**2)
        #print("homozygosity based theta :")
        thetaF = 0.5*(1/(F**2) - 1)
        print("thetaF : " + str(thetaF))
        Ne = thetaF/(4*Musat)
        lthetaF[i] = thetaF
        print("effective popualtion size : " + str(Ne))
        Vs = np.var(nuc[:, i+1])
        thetaV = 2*Vs
        print("thetaV : " +str(thetaV))
        Ne = thetaV/(4*Musat)
        print("effective popualtion size : " + str(Ne))
        lthetaV[i] = thetaV

    print("mean thetaF :", np.mean(lthetaF))
    print("mean thetaV :", np.mean(lthetaV))
    print("sex ratio :" +str(Nm/(Nm+Nf)))
    Ne = (Nm*Nf*4)/(Nm+Nf)
    print("population size estimated from sex ratio : " + str(Ne))

    #effective size estimation for cytoplasmic loci
    alleles, freqs  = np.unique(cyt, return_counts=True)
    freqs=freqs/N
    He = 1 - sum(freqs**2)
    thetaC = He/((1-He))
    print("cyt theta :" +str(thetaC))
    print("estimated effective pop size : " + str(thetaC/(4*Mucyt)))
    print("number of female:",str(Nf))
    print("estimated pop size:",str(Nf/2))

    f.close()
    return

try:
    mkdir("Simulations")
except:
    pass
Musat = 0.5e-2
N = 1000
print("theta theo : " + str(4*N*Musat))
run(N,[chr(i+65) for i in range(26)], Musat, 1e-3, 4500)

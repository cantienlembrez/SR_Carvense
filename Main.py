import numpy as np
from os import mkdir
from datetime import datetime

from BasicFunctions import *
from Models import *

def CorrThetaF(n, ThetaF):
    # biais correction see Xu and Fu 2004
    if ThetaF>15:
        A = (1.1675 + 3.3232/n + 63.698/(n**2))
        B = 0.2569
    else:
        A = (1.1313 + 3.4882/n + 28.2878/(n**2))
        B = 0.3998
    UnbiasedThetaF = ((np.sqrt(B**2 + 4*A*ThetaF) - B)/(2*A))**2
    return UnbiasedThetaF

def run(Model, Gmax, Dyn, N, IDmsats, Musat, Muiloc, Sm=0.5, em = 0.65, s = 0.44, d = 0.1, a = 5, g = 1.1):
    """
    Gmax: generation ending simulation
    Model: M1 Wright Fisher model, M2 male mortality, M3 Trioecy
    Dyn: a tupe (K, g) -> once Gmax reached saves the state of the simulation K times each g generations
    N: number of individuals,
    IDmsats: scaled position of microsatellites loci on chromosome (increasing order)
    Musat: microsatellites mutation rate
    Muiloc: mutation rate for infinite allele loci
    em: Proportion of male aborted pollen due to CMS
    Sm: Survial  rate of male model 2
    s: hermaphrodite selfing rate
    d: prop of selfing ovule aborted due to inbreeding depretion
    a: ratio  1  male pollen / 1 hermaphrodite pollen
    g: ratio 1 female ovules / 1 hermaphrodite ovules
    """

    N2 = N*2
    nMsats = len(IDmsats)

    #intialization
    if Model=="M1" or Model=="M2":
            nuc, cyt = init_m1m2(N, N2, nMsats)
            tmp_cyt = np.zeros(N,dtype=int)
    elif Model=="M3":
        HO = 1/g
        HP = 1/a
        MCMSP = 1-em
        nuc, cyt, Pollen, Ovules = init_CMS(N, N2, nMsats, HO, HP, MCMSP)
        tmp_cyt = np.zeros((N, 2),dtype=int)
    tmp_nuc = np.zeros((N2, nMsats+2),dtype=int)

    fpath = "Simulations/"+Model+"_"+datetime.now().strftime('%m_%d_%H:%M:%S')
    mkdir(fpath)
    f = open(fpath+"/parameter.txt", "a")
    f.write("Simulated generations:"+str(Gmax)+"\nParameters\nN:"+str(N)+"\nµ microsat:"+str(Musat)+"\nµ infinite allele loci:"+ str(Muiloc))
    if Model=="M2":
        f.write("\nSm:"+str(Sm))
    if Model=="M3":
        f.write("\ng:"+str(g)+"\na:"+str(a)+"\ns:"+str(s)+"\nd:"+str(d)+"\nem:"+str(em))
    f.write("\nInitial microsatellites state:"+str(nMsats)+"\n" +str(IDmsats)+"\n" + str(nuc[0,2:]))
    f.close()

    #generation loop
    if Model=="M1":
        for i in range(1, Gmax+1):
            gen_incr(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt)
            nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
        for dy in range(Dyn[0]*Dyn[1]):
            gen_incr(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt)
            nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
            if dy%Dyn[0]==0:
                np.save(fpath+"/nuc"+str(i+dy), nuc)
                np.save(fpath+"/cyt"+str(i+dy), cyt)
    elif Model=="M2":
        for i in range(1, Gmax+1):
            gen_incr_m2(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt, Sm)
            nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
        for dy in range(Dyn[0]*Dyn[1]):
            gen_incr_m2(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt, Sm)
            nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
            if dy%Dyn[0]==0:
                np.save(fpath+"/nuc"+str(i+dy), nuc)
                np.save(fpath+"/cyt"+str(i+dy), cyt)
    elif Model=="M3":
        for i in range(1, Gmax+1):
            Ovules, Pollen = gen_incr_CMS(N, N2, nMsats, Musat, Muiloc, HO, HP, MCMSP, s, d, nuc, cyt, Ovules, Pollen, tmp_nuc, tmp_cyt)
            nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
        for dy in range(Dyn[0]*Dyn[1]):
            Ovules, Pollen = gen_incr_CMS(N, N2, nMsats, Musat, Muiloc, HO, HP, MCMSP, s, d, nuc, cyt, Ovules, Pollen, tmp_nuc, tmp_cyt)
            nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
            if dy%Dyn[0]==0:
                np.save(fpath+"/nuc"+str(i+dy), nuc)
                np.save(fpath+"/cyt"+str(i+dy), cyt)
    return

try:
    mkdir("Simulations")
except:
    pass

Mu = 1e-3
N = 400
run("M3", 1600, (10, 10), N,[chr(i+65) for i in range(1)], Mu, Mu)

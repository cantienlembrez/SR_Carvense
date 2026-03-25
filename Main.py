import numpy as np
from os import mkdir

from BasicFunctions import *
from Models import *
import sys

F_ID = 0

def run(Model, Gmax, Dyn, N, IDmsats, Musat, Muiloc, Sm=None, em=None, s=None, d=None, a=None, g=None):
    """
    Gmax: generation ending simulation
    Model: M1 Wright Fisher model, M2 male mortality, M3 Trioecy
    Dyn: a tuple (g, K) -> once Gmax reached saves the state of the simulation K times each g generations
    N: number of individuals,
    IDmsats: ID of microsatellites loci on chromosome (increasing order)
    Musat: microsatellites mutation rate
    Muiloc: mutation rate for infinite allele loci
    em: Proportion of male aborted pollen due to CMS
    Sm: Survial  rate of male model 2
    s: hermaphrodite selfing rate
    d: prop of selfing ovule aborted due to inbreeding depretion
    a: ratio  1  male pollen / 1 hermaphrodite pollen
    g: ratio 1 female ovules / 1 hermaphrodite ovules
    """
    global F_ID

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

    fpath = "Simulations/"+Model+"_"+"N"+str(N)+"_"+str(F_ID)
    F_ID+=1
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
            gen_incr_m1(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt)
            nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
        for dy in range(Dyn[0]*Dyn[1]):
            gen_incr_m1(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt)
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
nMsats = 20


MODEL = sys.argv[1]
N = int(sys.argv[2])
NbSave = int(sys.argv[3])
Int = int(sys.argv[4])
Nbreplicates = int(sys.argv[5])
SEED = int(sys.argv[6])

Sm, em, s, d, a, g = None, None, None, None, None, None

for A in range(6, len(sys.argv)):
    Arg = sys.argv[A].split(":")
    if Arg[0]=="Sm": Sm = float(Arg[1])
    elif Arg[0]=="em": em = float(Arg[1])
    elif Arg[0]=="s": s = float(Arg[1])
    elif Arg[0]=="d": d = float(Arg[1])
    elif Arg[0]=="a": a = float(Arg[1])
    elif Arg[0]=="g": g = float(Arg[1])


np.random.seed(SEED)
for replicate in range(Nbreplicates):
    run(MODEL, 6*N, (Int, NbSave), N,[chr(i+65) for i in range(nMsats)], Mu, Mu, Sm, em, s, d, a, g)

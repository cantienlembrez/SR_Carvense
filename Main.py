import numpy as np
from os import mkdir
import sys

from BasicFunctions import *
from Models import *

F_ID = 0

def run(Model, CP, Gmax, Dyn, N, IDmsats, Musat, Muiloc, Sm=None, em=None, s=None, d=None, a=None, g=None, MaxAge = None, p = None, c = None, K = None):
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
    CP: True or false, True -> perennity/clonality simulated
    """

    global F_ID

    N2 = N*2
    nMsats = len(IDmsats)
    mut_ID = [1, 1] # ID for infinite allele loci [ Nuc, Cyt ]

    #intialization
    if Model=="M1" or Model=="M2":
        nuc, cyt = init_m1m2(N, N2, nMsats)
        tmp_cyt = np.zeros(N,dtype=int)
        Pollen, Ovules = None, None
        parameters = []
        sex_func = sex_m1
        if Model=="M2":
            parameters = [Sm]
            sex_func = sex_m2
    elif Model=="M3":
        HO = 1/g
        HP = 1/a
        MCMSP = 1-em
        nuc, cyt, Pollen, Ovules = init_CMS(N, N2, nMsats, HO, HP, MCMSP)
        tmp_cyt = np.zeros((N, 2),dtype=int)
        parameters = [HO, HP, MCMSP, s, d]
        sex_func = sex_CMS
    tmp_nuc = np.zeros((N2, nMsats+2),dtype=int)
    print("Init", F_ID)
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

    cp_func = ClonalityPerenity
    if CP==1:
        f.write("\nCP \np:"+str(p)+"\nc:"+str(c))
        Age = np.zeros(N, dtype=int)
        cp_param = [MaxAge, p, c]
        cl_func = ClonalityUB
    elif CP==2:
        f.write("\nCPB \np:"+str(p)+"\nc:"+str(c)+"\nK:"+str(K))
        Age = np.zeros(N, dtype=int)
        cp_param = [MaxAge, p, c, c*K]
        cl_func = ClonalityB
    else:
        f.write("\nNo CP")
        Age = None
        cp_param = None
        cp_func = No_CP
        cl_func = None
    f.close()

    #generation loop
    print("RUN")
    for i in range(1, Gmax+1):
        Nsurv = cp_func(cp_param, N, nuc, cyt, tmp_nuc, tmp_cyt, Age, cl_func)
        Ovules, Pollen = sex_func(parameters, N, N2, Nsurv, nMsats, Musat, Muiloc, mut_ID, nuc, cyt, tmp_nuc, tmp_cyt, Ovules, Pollen)
        nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
        # test male extinction
        if np.sum(nuc[1:N2:2,0])==0:
            f = open(fpath+"/m_extinction"+str(i)+".txt", "a")
            f.close()
            return
    for dy in range(Dyn[0]*Dyn[1]): #loop saving the system state after Gmax generations
        Nsurv = cp_func(cp_param, N, nuc, cyt, tmp_nuc, tmp_cyt, Age, cl_func)
        Ovules, Pollen = sex_func(parameters, N, N2, Nsurv, nMsats, Musat, Muiloc, mut_ID, nuc, cyt, tmp_nuc, tmp_cyt, Ovules, Pollen)
        nuc, cyt = np.copy(tmp_nuc), np.copy(tmp_cyt)
        # test male extinction
        if np.sum(nuc[1:N2:2,0])==0:
            f = open(fpath+"/m_extinction"+str(i)+".txt", "a")
            f.close()
        if dy%Dyn[0]==0:
            np.save(fpath+"/nuc"+str(i+dy), nuc)
            np.save(fpath+"/cyt"+str(i+dy), cyt)
            if CP>0:
                np.save(fpath+"/Age"+str(i+dy), Age)
    return

try:
    mkdir("Simulations")
except:
    pass

Mu = 1e-3
nMsats = 10

MODEL = sys.argv[1]
CP = int(sys.argv[2])
N = int(sys.argv[3])
NbSave = int(sys.argv[4])
Int = int(sys.argv[5])
Nbreplicates = int(sys.argv[6])
SEED = int(sys.argv[7])

Sm, em, s, d, a, g, MaxAge, p, c, K = None, None, None, None, None, None, None, None, None, None
Gmax = 6*N

for A in range(7, len(sys.argv)):
    Arg = sys.argv[A].split(":")
    if Arg[0]=="Sm": Sm = float(Arg[1])
    elif Arg[0]=="em": em = float(Arg[1])
    elif Arg[0]=="s": s = float(Arg[1])
    elif Arg[0]=="d": d = float(Arg[1])
    elif Arg[0]=="a": a = float(Arg[1])
    elif Arg[0]=="g": g = float(Arg[1])
    elif Arg[0]=="p": p = float(Arg[1])
    elif Arg[0]=="c": c = float(Arg[1])
    elif Arg[0]=="MaxAge": MaxAge = float(Arg[1])
    elif Arg[0]=="K": K = float(Arg[1])
    elif Arg[0]=="Gmax": Gmax = int(Arg[1])

np.random.seed(SEED)
for replicate in range(Nbreplicates):
    run(MODEL, CP, Gmax, (Int, NbSave), N,[chr(i+65) for i in range(nMsats)], Mu, Mu, Sm, em, s, d, a, g, MaxAge, p, c, K)

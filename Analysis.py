import numpy as np

def CorrThetaF(n, ThetaF):
    # biais correction for thetaF see Xu and Fu 2004
    if ThetaF>15:
        A = (1.1675 + 3.3232/n + 63.698/(n**2))
        B = 0.2569
    else:
        A = (1.1313 + 3.4882/n + 28.2878/(n**2))
        B = 0.3998
    UnbiasedThetaF = ((np.sqrt(B**2 + 4*A*ThetaF) - B)/(2*A))**2
    return UnbiasedThetaF


def EffectiveSize(dirname, nucfile, cytfile):
    Model = dirname[0:2]
    dirname = "Simulations/"+dirname
    fp = open(dirname+"parameter.txt", "r")
    for l in fp:
        splitline = l.split(":")
        if splitline[0]=="N":
            N=int(splitline[1][0:-1])
        elif splitline[0]=="µ microsat":
            Musat=float(splitline[1][0:-1])
        elif splitline[0]=="µ infinite allele loci":
            Muiloc=float(splitline[1][0:-1])
        elif splitline[0]=="s":
            s=float(splitline[1][0:-1])
        elif splitline[0]=="Initial microsatellites state":
            nMsats=int(splitline[1][0:-1])
    nuc = np.load(dirname+nucfile)
    cyt = np.load(dirname+cytfile)
    print("nuclear Effective size :")
    if Model=="M1" or Model=="M2":
        Nm = np.sum(nuc[:, 0])
        Nf = N - Nm
        print("Sex Ratio :", Nm/(N))
        print("Effective Pop size :", 4*Nm*Nf/N)
    elif Model=="M3":
        Nh, Nm, Nf, Nmcms = 0, 0, 0, 0
        for i in range(N):
            if nuc[i*2+1,0]==1:
                if cyt[i,0]==1:
                    Nmcms+=1
                Nm+=1
            if nuc[i*2+1,0]==0 and cyt[i,0]==1:
                Nf+=1
            if nuc[i*2+1,0]==0 and cyt[i,0]==0:
                Nh+=1
        print("Number of males", Nm, "("+str(Nmcms)+" with CMS)")
        print("Number of hermaphrodites", Nh)
        print("Number of female", Nf)
        print("Effective Pop size :", 4*(Nm + (1-s)*Nh)*(Nf + (1-s)*Nh)/(Nm+Nf+2*(1-s)*Nh)) #warning
    # nuclear infinite allele locus
    freqs = np.unique(nuc[:, 1], return_counts=True)[1]/(N*2)
    He = 1 - np.sum(freqs**2)
    theta = He/(1-He)
    Ne = theta/(4*Muiloc)
    print("estimation with infinite allele locus", Ne)
    NeF = [None] * nMsats
    NeV = [None] * nMsats
    for i in range(nMsats):
        freqs = np.unique(nuc[:, i+2], return_counts=True)[1]/(N*2)
        F = np.sum(freqs**2)
        thetaF = 0.5*(1/(F**2) - 1)
        NeF[i] = CorrThetaF(N, thetaF)/(4*Musat)
        thetaV = np.var(nuc[:, i+2])*2
        NeV[i] = thetaV/(4*Musat)
    print("Mean estimation with microsat thetaF", np.mean(NeF))
    print("Mean estimation with microsat thetaV", np.mean(NeV))
    print("--------------------------------------------\ncytoplasmic effective size :")
    if Model=="M1" or Model=="M2":
        print("estimation with number of female", Nf/2)
        freqs = np.unique(cyt, return_counts=True)[1]/(N*2)
    if Model=="M3":
        #todo
        freqs = np.unique(cyt[:,1], return_counts=True)[1]/(N*2)
    He = 1 - np.sum(freqs**2)
    theta = He/(1-He)
    Ne = theta/(2*Muiloc)
    print("estimation with infinite allele locus", Ne)
    return


EffectiveSize("M2_03_21_11:00:49/", "nuc490.npy", "cyt490.npy")

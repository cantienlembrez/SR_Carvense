import numpy as np
import matplotlib.pyplot as plt

from AnalysisFuncs import *

path = "[...]"

N=1000
He_null6Ncp = np.zeros(100, float)
He_null7Ncp = np.zeros(100, float)
He_null8Ncp = np.zeros(100, float)
Nall_nullcp = np.zeros(100, float)
Ne_nullcp = np.zeros((100, 3), float)
S_nullcp = np.zeros(100, dtype=float)
Gen_nullcp = np.zeros((100, 3), float)
List_AllFreqs_nullcp = []
for i in range(100):
    dirname = path+"/M1/M1_N"+str(N)+"_"+str(i)+"/"
    fp = open(dirname+"parameter.txt", "r")
    for l in fp:
        splitline = l.split(":")
        if splitline[0]=="N":
            N=int(splitline[1][0:-1])
        elif splitline[0]=="µ microsat":
            Musat=float(splitline[1][0:-1])
        elif splitline[0]=="µ infinite allele loci":
            Muiloc=float(splitline[1][0:-1])
        elif splitline[0]=="Initial microsatellites state":
            nMsats=int(splitline[1][0:-1])
    nuc = np.load(dirname+"nuc"+str(6*N)+".npy")
    He_null6Ncp[i] = (1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
    nuc = np.load(dirname+"nuc"+str(6*N+1000)+".npy")
    He_null7Ncp[i] = (1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
    nuc = np.load(dirname+"nuc"+str(6*N+2000)+".npy")
    He_null8Ncp[i] = (1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
    Ne_nullcp[i,:] = NucEffectiveSize(N, "M1", nuc, Musat, Muiloc, nMsats)
    S_nullcp[i] = np.sum(nuc[:,0])/N
    Gen_nullcp[i, : ] = Unique_genotype(N, nuc)
    Nall_nullcp[i] = NAll(nuc)
    for k in range(nMsats):
        List_AllFreqs_nullcp.append(np.unique(nuc[:,k+2], return_counts = True)[1]/(N*2))

H = np.mean(He_null6Ncp)
print("Ne null estimation with inf allele locus", H/((1-H)*4*Muiloc))

#### Figures ####

#Equilibre
plt.figure("equilibre")
ts = np.array([np.mean(He_null6Ncp),np.mean(He_null7Ncp), np.mean(He_null8Ncp)])
stdts = np.array([np.std(He_null6Ncp),np.std(He_null7Ncp), np.std(He_null8Ncp)])
plt.errorbar(np.array([6000, 7000, 8000]), ts, stdts, marker='^')

plt.show()

# Ne
plt.figure("thetaF")
counts, bins = np.histogram(np.concat((Ne_nullcp[:,1], Ne_null[:,1])))
plt.hist(Ne_nullcp[:,1], bins = bins, label="c=0.1 p=1 (100 replicates)", color ="red")
plt.hist(Ne_null[:,1], bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Ne estimated from microsat homozigoty")
plt.legend()
plt.show()

# Ne
plt.figure("thetaF")
counts, bins = np.histogram(np.concat((Ne_nullcp[:,1], Ne_null[:,1])))
plt.hist(Ne_nullcp[:,1], bins = bins, label="c=0.1 p=1 (100 replicates)", color ="red")
plt.hist(Ne_null[:,1], bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Ne estimated from microsat homozigoty")
plt.legend()
plt.show()

plt.figure("thetaV")
counts, bins = np.histogram(np.concat((Ne_nullcp[:,2], Ne_null[:,2])))
plt.hist(Ne_nullcp[:,2], bins = bins, label="c=0.1 p=1 (100 replicates)", color ="red")
plt.hist(Ne_null[:,2], bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Ne estimated from microsat variance in size")
plt.legend()
plt.show()

# Nb all
plt.figure("NbAll")
counts, bins = np.histogram(np.concat((Nall_nullcp, Nall_null)))
plt.hist(Nall_nullcp, bins = bins, label="c=0.1 p=1 (100 replicates)", color ="red")
plt.hist(Nall_null, bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Mean allele number per loci")
plt.legend()
plt.show()

#SFS
List_AllFreqs_nullcp = np.concat(List_AllFreqs_nullcp)
List_AllFreqs_null = np.concat(List_AllFreqs_null)

list_global = np.concat((List_AllFreqs_nullcp, List_AllFreqs_null))

plt.figure("SFS")
counts, bins = np.histogram(list_global)
plt.hist(List_AllFreqs_nullcp, bins = bins, label="c=0.1 p=1 (100 replicates)", color ="red")
plt.hist(List_AllFreqs_null, bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("counts of alleles frequencies")
plt.legend()
plt.show()

plt.figure("SFS 2")
counts, bins = np.histogram(list_global[list_global<=0.05])
plt.hist(List_AllFreqs_nullcp, bins = bins, label="c=0.1 p=1 (100 replicates)", color ="red")
plt.hist(List_AllFreqs_null, bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("counts of alleles frequencies")
plt.legend()
plt.show()

# Genotypes
plt.figure("Nb genotypes")
plt.title("nb of unique genotypes")
counts, bins = np.histogram(np.concat((Gen_null[:,0], Gen_nullcp[:,0])))
plt.hist(Gen_nullcp[:,0], bins=bins, label="c=0.1 p=1 (100 replicates)", color ="red")
plt.hist(Gen_null[:,0], bins=bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.legend()
plt.show()

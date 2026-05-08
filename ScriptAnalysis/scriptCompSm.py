import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *

from AnalysisFuncs import *

path_paramset = "[...]"
#path to folder_paramset/
# M3/
#   All M3 simulations
# M1eq/
#   All equivalent M1 simulations
Npair = (5000, 49)
Musat = 1e-3
Muiloc =1e-3
nMsats = 10


He_mrt = []
S_mrt = []
Hec_mrt = []
Ne_mrt = []
Nall_mrt = []
Nall_infmrt = []
List_AllFreqs_mrt = []
Fis_mrt = []
Gen_mrt = np.zeros((100,3), dtype=float)
N = Npair[0]
for i in range(100):
    dirname = path_paramset+ "M2/M2_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)

    He_mrt.append(1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
    Hec_mrt.append(1 - np.sum((np.unique(cyt, return_counts = True)[1]/(N))**2))
    Ne_mrt.append(NucEffectiveSize(N, "M2", nuc, Musat, Muiloc, nMsats))
    Nall_infmrt.append(len(np.unique(nuc[:,1])))
    Nall_mrt.append(NAll(nuc))
    for k in range(nMsats):
        List_AllFreqs_mrt.append(np.unique(nuc[:,k+2], return_counts = True)[1]/(N*2))
    Fis_mrt.append(Fis(N, nuc))
    Gen_mrt[i, :] = Unique_genotype(N, nuc)
    S_mrt.append(SexNumber("M1", N, nuc))

SR_mean = np.mean(S_mrt, axis=0)
print("Nm, Nf :", SR_mean/Npair[0])
#expected = np.array((1/3, 2/3))
expected = np.array((3/7, 4/7))
expected = expected * Npair[0]
X2 = chisquare(SR_mean, expected)
print(X2)

He_null = []
Hec_null = []
Ne_null = []
Nall_null = []
Nall_infnull = []
List_AllFreqs_null = []
Fis_null = []
Gen_null = np.zeros((100,3), dtype=float)
N = Npair[1]
for i in range(100):
    dirname = path_paramset+"M1eq/M1_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)


    He_null.append(1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
    Hec_null.append(1 - np.sum((np.unique(cyt, return_counts = True)[1]/(N))**2))
    Ne_null.append(NucEffectiveSize(N, "M1", nuc, Musat, Muiloc, nMsats))
    Nall_infnull.append(len(np.unique(nuc[:,1])))
    Nall_null.append(NAll(nuc))
    for k in range(nMsats):
        List_AllFreqs_null.append(np.unique(nuc[:,k+2], return_counts = True)[1]/(N*2))
    Fis_null.append(Fis(N, nuc))
    Gen_null[i, :] = Unique_genotype(N, nuc)

Ne_mrt = np.array(Ne_mrt)
Ne_null = np.array(Ne_null)

# Stats Ne
print("-------------------------Ne-------------------------------")
print("thetaF Ne mrt : ", np.mean(Ne_mrt[:,0]), np.std(Ne_mrt[:,0]))
#print("estimators Ne null : ", np.mean(Ne_null, axis=0))
print(mannwhitneyu(Ne_mrt[:,0], Ne_null[:,0]))

print("thetaV Ne mrt : ", np.mean(Ne_mrt[:,1]), np.std(Ne_mrt[:,1]))
print(mannwhitneyu(Ne_mrt[:,1], Ne_null[:,1]))

print("He_mrt :", np.mean(He_mrt), np.std(He_null))
H = np.mean(He_mrt)
print("Ne mrt estimation with inf allele locus", H/((1-H)*4*Muiloc))
print(mannwhitneyu(He_mrt, He_null))

print("For cytoplamsic locus : ")
H = np.mean(Hec_mrt)
print("mrt mean H", H, "std : ", np.std(Hec_mrt), "eff size", H/((1-H)*4*Muiloc))
H = np.mean(Hec_null)
print("null mean H", H, "std : ", np.std(Hec_null), "eff size", H/((1-H)*4*Muiloc))

# Stats Fis
print("-------------------------Fis------------------------------")
print("Fis mrt : ", np.mean(Fis_mrt), "Fis_null :",np.mean(Fis_null))
print(mannwhitneyu(Fis_mrt, Fis_null))
print("Fis avec 0 : ", mannwhitneyu(Fis_mrt, 0))
# Stats Nball
print("-------------------------Nall------------------------------")
print("Nall mrt : ", np.mean(Nall_mrt), np.std(Nall_mrt), "Nall null :",np.mean(Nall_null), np.std(Nall_null))
print(mannwhitneyu(Nall_mrt, Nall_null))


print("inf all locus")
print("Nall mrt : ", np.mean(Nall_infmrt), np.std(Nall_infmrt), "Nall null :",np.mean(Nall_infnull), np.std(Nall_infnull))
print("difference between mean at inf loc", np.mean(Nall_infmrt) -  np.mean(Nall_infnull))
print(mannwhitneyu(Nall_infmrt, Nall_infnull))

#### Figures ####

#thetaF
plt.figure("thetaF")
bins = np.histogram_bin_edges(np.concat((Ne_mrt[:,1], Ne_null[:,1])), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_null = np.histogram(Ne_null[:,1], bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

counts_mrt = np.histogram(Ne_mrt[:,1], bins=bins)[0]
counts_mrt = counts_mrt/np.sum(counts_mrt)

plt.bar(x, counts_mrt, w, alpha=0.7, label="Diff mort model (100 replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Ne estimated from microsat homozigoty")
plt.legend()
plt.show()

#thetaV
plt.figure("thetaV")
bins = np.histogram_bin_edges(np.concat((Ne_mrt[:,2], Ne_null[:,2])), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_null = np.histogram(Ne_null[:,2], bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

counts_mrt = np.histogram(Ne_mrt[:,2], bins=bins)[0]
counts_mrt = counts_mrt/np.sum(counts_mrt)

plt.bar(x, counts_mrt, w, alpha=0.7, label="Diff mort model (100 replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Ne estimated from microsat variance in size")
plt.legend()
plt.show()


# Nb all
plt.figure("NbAll")
bins = np.histogram_bin_edges(np.concat((Nall_null, Nall_mrt)), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_mrt = np.histogram(Nall_mrt, bins=bins)[0]
counts_mrt = counts_mrt/np.sum(counts_mrt)

counts_null = np.histogram(Nall_null, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

plt.bar(x, counts_mrt, w, alpha=0.7, label="Diff mort model (100 replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Mean number of allele")
plt.legend()
plt.show()


plt.figure("NbAll")
bins = np.histogram_bin_edges(np.concat((Nall_infnull, Nall_infmrt)), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_mrt = np.histogram(Nall_infmrt, bins=bins)[0]
counts_mrt = counts_mrt/np.sum(counts_mrt)

counts_null = np.histogram(Nall_infnull, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

plt.bar(x, counts_mrt, w, alpha=0.7, label="Diff mort model (100 replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Mean number of allele")
plt.legend()
plt.show()

# SFS
List_AllFreqs_mrt = np.concat(List_AllFreqs_mrt)
List_AllFreqs_null = np.concat(List_AllFreqs_null)

list_global = np.concat((List_AllFreqs_mrt, List_AllFreqs_null))

plt.figure("SFS")
bins = np.histogram_bin_edges(list_global, bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]


counts_mrt = np.histogram(List_AllFreqs_mrt, bins=bins)[0]
counts_mrt = counts_mrt/np.sum(counts_mrt)

counts_null = np.histogram(List_AllFreqs_null, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

plt.bar(x, counts_mrt, w, alpha=0.7, label="Diff mort model (100 replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("counts of alleles frequencies")
plt.legend()
plt.show()

plt.figure("SFS 2")
bins = np.histogram_bin_edges(list_global[list_global<=0.1], bins=8)
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]


counts_mrt = np.histogram(List_AllFreqs_mrt, bins=bins)[0]
counts_mrt = counts_mrt/np.sum(counts_mrt)

counts_null = np.histogram(List_AllFreqs_null, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

plt.bar(x, counts_mrt, w, alpha=0.7, label="Diff mort model (100 replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("counts of alleles frequencies")
plt.legend()
plt.show()

# Fis
plt.figure("Fis")
bins = np.histogram_bin_edges(np.concat((Fis_null, Fis_mrt)), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_mrt = np.histogram(Fis_mrt, bins=bins)[0]
counts_mrt = counts_mrt/np.sum(counts_mrt)

counts_null = np.histogram(Fis_null, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

plt.bar(x, counts_mrt, w, alpha=0.7, label="Diff mort model (100 replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Fis")
plt.legend()
plt.show()


#### Genotypes ####

plt.figure("Nb genotypes")
plt.title("nb of unique genotypes")
counts, bins = np.histogram(np.concat((Gen_null[:,0], Gen_mrt[:,0])))
plt.hist(Gen_mrt[:,0], label="Diff mort model (100 replicates)", color ="red")
plt.hist(Gen_null[:,0], alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.legend()
plt.show()

plt.figure("proportion unique genotypes")
plt.title("proportion of unique genotypes")
counts, bins = np.histogram(np.concat((Gen_null[:,1], Gen_mrt[:,1])))
plt.hist(Gen_mrt[:,1],  bins = bins, label="Diff mort model (100 replicates)", color ="red")
plt.hist(Gen_null[:,1], bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.legend()
plt.show()

plt.figure("Clonal eveness")
plt.title("clonal eveness")
counts, bins = np.histogram(np.concat((Gen_null[:,2], Gen_mrt[:,2])))
plt.hist(Gen_mrt[:,2],  bins = bins, label="Diff mort model (100 replicates)", color ="red")
plt.hist(Gen_null[:,2], bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.legend()
plt.show()

plt.figure("prop unique Genotypes function of mean number of alleles")
plt.plot(Nall_mrt/np.max(Nall_mrt), Gen_mrt[:,1], "o", color ="red")
plt.plot(Nall_null/np.max(Nall_null), Gen_null[:,1], "o", color ="grey")
plt.xlabel("Normalized mean allele number per loci")
plt.ylabel("proportion of unique")
plt.show()


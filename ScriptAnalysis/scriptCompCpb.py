import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *

from AnalysisFuncs import *

path = "[...]"
#path to the dir containing all simulations dirs

N=1000
Musat = 1e-3
Muiloc = 1e-3
nMsats = 10
Ne_cpb = np.zeros((100, 2), float)
S_cpb = np.zeros((100, 2), dtype=float)
Gen_cpb = np.zeros((100, 3), float)
Eveness = np.zeros((100, 6), float)
He_cpb = np.zeros((100), float)
Fis_cpb = np.zeros((100), float)
Nall_cpb = np.zeros(100, float)
He_cyt = np.zeros(100, float)
for i in range(100):
    dirname = path+"/M1_N"+str(N)+"_"+str(i)+"/"
    nuc = np.load(dirname+"nuc"+str(6*N)+".npy")
    cyt = np.load(dirname+"cyt"+str(6*N)+".npy")
    He_cpb[i] = (1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
    Ne_cpb[i,:] = NucEffectiveSize(N,"M1", nuc, Musat, Muiloc, nMsats)
    Nall_cpb[i] = NAll(nuc)
    Fis_cpb[i] = Fis(N, nuc)
    S_cpb[i,] = np.array([np.sum(nuc[:,0]), N-np.sum(nuc[:,0])])/N
    Gen_cpb[i] = Unique_genotype(N, nuc)
    for nbmarker in [5, 6, 7, 8, 9, 10]:
        Eveness[i, 10-nbmarker] = Unique_genotype(N, nuc, loci_to_compare = range(2,nbmarker+2))[2]
    He_cyt[i] = (1 - np.sum((np.unique(cyt, return_counts = True)[1]/(N))**2))
print("median Fis != 0", mannwhitneyu(Fis_cpb, 0))

print("thetaF Ne cpb : ", np.mean(Ne_cpb[:,0]), np.std(Ne_cpb[:,0]))


print("thetaV Ne cpb : ", np.mean(Ne_cpb[:,1]), np.std(Ne_cpb[:,1]))


print("He_cpb :", np.mean(He_cpb), np.std(He_cpb))
H = np.mean(He_cpb)
print("Ne cpb estimation with inf allele locus", H/((1-H)*4*Muiloc))


print("Nm Nf : ", np.mean(S_cpb, axis=0))
H = np.mean(He_cyt)
print("clon mean H", H, "std : ", np.std(He_cyt), "eff size", H/((1-H)*4*Muiloc))


#### Génotypes ####
Gen_cp = np.zeros((100, 3), float)
Nall_cp = np.zeros(100, float)
for i in range(100):
    dirname = "/home/clembrez/Documents/Master/M1/Stage/Donnees_Simulees/TextClonalityMaxAge2C0.1P1/M1/M1_N"+str(N)+"_"+str(i)+"/"
    nuc = np.load(dirname+"nuc"+str(6*N)+".npy")
    Gen_cp[i] = Unique_genotype(N, nuc)
    Nall_cp[i] = NAll(nuc)



#### Effet Nb alleles ####
mod1 = linregress(Nall_cpb, Gen_cpb[:,0])
pred = mod1.slope * Nall_cpb + mod1.intercept
residuals = Gen_cpb[:,0] - pred

# plt.figure("homoscedastitité")
# plt.scatter(pred, residuals)
# plt.axhline(0, color='red')
# plt.show()

shapiro(residuals)
# A peu pres OK
print(mod1.pvalue)

mod2 = linregress(Nall_cp, Gen_cp[:,0])
pred = mod2.slope * Nall_cp + mod2.intercept
residuals = Gen_cp[:,0] - pred

# plt.figure("homoscedastitité2")
# plt.scatter(pred, residuals)
# plt.axhline(0, color='red')
# plt.show()

shapiro(residuals)
# pareil A peu pres OK
print(mod2.pvalue)

plt.rcParams.update({'font.size': 18})
plt.figure("Nb Genotype function of mean number of alleles")
plt.subplots_adjust(left=0.1, right=0.90, top=1, bottom=0.15)
plt.plot(Nall_cpb/np.max(Nall_cpb), Gen_cpb[:,0], "o", color ="red", label = "Simulations avec biais de clonalité (pente : $0.22$, $t=8.00$, 98 ddl, $p=2.5\cdot 10^{-12}$)")
plt.plot(Nall_cp/np.max(Nall_cp), Gen_cp[:,0], "o", color ="grey", label = "Simulations sans biais de clonalité (pente : $0.18$, $t=8.04$, 98 ddl, $p=2.0\cdot 10^{-12}$)")

plt.plot(np.linspace(np.min(Nall_cpb), np.max(Nall_cpb), 100)/np.max(Nall_cpb), mod1.slope * np.linspace(np.min(Nall_cpb), np.max(Nall_cpb), 100) + mod1.intercept, color="red")
plt.plot(np.linspace(np.min(Nall_cp), np.max(Nall_cp), 100)/np.max(Nall_cp), mod2.slope * np.linspace(np.min(Nall_cp), np.max(Nall_cp), 100) + mod2.intercept, color="grey")

plt.ylabel("Proportion de génotypes uniques")
plt.xlabel("Nombre d'allèles par locus normalisé")
plt.legend()
plt.show()




### Distributions ###
## Prop G uniques
bins = np.histogram_bin_edges(np.concatenate((Gen_cp[:,0], Gen_cpb[:,0])), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]


counts_cp = np.histogram(Gen_cp[:,0], bins=bins)[0]
counts_cp = counts_cp/np.sum(counts_cp)

counts_cpb = np.histogram(Gen_cpb[:,0], bins=bins)[0]
counts_cpb = counts_cpb/np.sum(counts_cpb)

plt.bar(x, counts_cpb, w, alpha=0.7, label="Générations chevauchantes avec biais de clonalité", color ="red")
plt.bar(x, counts_cp, w,  alpha=0.7, label="Générations chevauchantes sans biais", color="grey")
plt.legend()
plt.show()

print(mannwhitneyu(Gen_cp[:,0],Gen_cpb[:,0]))

## eveness

pearsonr(Eveness, np.repeat(Nall_cpb, 6, axis=0).reshape(100, 6), alternative='two-sided', method=None, axis=0) # when it comes under 9 loci genotypes are limited by allele number

plt.rcParams.update({'font.size': 18})

plt.subplots_adjust(left=0.05, right=0.95, top=1, bottom=0.15)
bins = np.linspace(0.7, 1, 10)
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]
b_w = (w*0.9)/5

Nb = [10, 9, 8, 7, 6]
for i in range(Eveness.shape[1] - 1):
    counts_cpb = np.histogram(Eveness[:,i], bins=bins)[0]
    counts_cpb = counts_cpb/np.sum(counts_cpb)
    plt.bar(x + b_w*(i-2), counts_cpb, b_w, label=str(Nb[i])+" Marqueurs")

bins[0] = 0
labels = [f"{bins[i]:.2f}-{bins[i+1]:.2f}" for i in range(len(bins)-1)]
plt.xticks(x, labels, rotation=45)
plt.xlabel("$E_{1/D}$")
plt.ylabel("Frequences")
plt.legend()
plt.show()

print(mannwhitneyu(Gen_cp[:,2],Gen_cpb[:,2]))
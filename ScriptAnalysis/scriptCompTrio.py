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
Npair = (1000, 132)
Musat = 1e-3
Muiloc =1e-3
nMsats = 10

Dist_mf_trio = []
He_trio = []
Hall = []
Ne_trio = []
Nall_trio = []
Nall_inftrio = []
List_AllFreqs_trio = []
Fis_trio = []
Gen_trio = []
S_trio = []
#Nball_sex = []
Fis_sex = []
Dist_mat = []
Dist_Distrib = []
Hec_trio = []
non_fix = 0
N = Npair[0]
for i in range(100):
    dirname = path_paramset+"M3/M3_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)

    #loading 6N state
    nuc = np.load(dirname+"nuc"+str(6*N)+".npy")
    cyt = np.load(dirname+"cyt"+str(6*N)+".npy")

    S_trio.append(SexNumber("M3", N, nuc, cyt = cyt))
    if S_trio[i][0] != 0: #dioecy not fixed
        non_fix+=1
        # Nuclear stuff
        He_trio.append(1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
        Hall.append(heterozigoty(N, nuc))
        Ne_trio.append(NucEffectiveSize(N, "M3", nuc, Musat, Muiloc, nMsats))
        Nall_trio.append(NAll(nuc))
        Nall_inftrio.append(len(np.unique(nuc[:,1])))
        for k in range(nMsats):
            List_AllFreqs_trio.append(np.unique(nuc[:,k+2], return_counts = True)[1]/(N*2))
        Fis_trio.append(Fis(N, nuc))
        Gen_trio.append(Unique_genotype(N, nuc))

        #Cytoplasmic stuff
        w_cms = cyt[cyt[:,0]==0,]
        c_cms = cyt[cyt[:,0]==1,]
        Hec_trio.append(1 - np.sum((np.unique(cyt[:,1], return_counts = True)[1]/(N))**2))
        #Per sex class comparaison
        #Nball_sex.append(NAll(nuc, True, [True, True, cyt]))
        h, m, c, f = split_nuc(nuc, True, cyt)
        if S_trio[i][1]!=0:
            Fis_sex.append([
            Fis(S_trio[i][0], h),
            Fis(S_trio[i][1], m),
            Fis(S_trio[i][2], c),
            Fis(S_trio[i][3], f)])
            #Dist, Distrib =  Dist_sexmatrix(nuc, cyt, True, 1000)
            #Dist_mat.append(Dist)
            #Dist_Distrib.append(Distrib)
            Dist_mf_trio.append(Dist_mf(nuc, True, cyt))
Dist_mf_trio = np.array(Dist_mf_trio)

S_trio = np.array(S_trio)

SR_mean = np.mean(S_trio[S_trio[:,0]!=0], axis=0)
print("Nh, Nm, Nc, Nf :", SR_mean/N)
#expected = np.array([0.073091,   0.01686715, 0.34126569, 0.56877615])
expected = np.array([0.06564365, 0.00724638, 0.35024155, 0.57686843])
expected = expected * Npair[0]
print("expected       :", expected)
X2 = chisquare(SR_mean, expected)
print(X2)

Dist_mf_null = []
He_null = []
Ne_null = []
Nall_null = []
Nall_infnull = []
List_AllFreqs_null = []
Fis_null = []
Gen_null = []
Hec_null = []
N = Npair[1]
for i in range(100):
    dirname = path_paramset + "M1eq/M1_N"+str(N)+"_"+str(i)+"/"
    nuc, cyt = load_sim(dirname, N)

    He_null.append(1 - np.sum((np.unique(nuc[:,1], return_counts = True)[1]/(N*2))**2))
    Ne_null.append(NucEffectiveSize(N, "M1", nuc, Musat, Muiloc, nMsats))
    Hec_null.append(1 - np.sum((np.unique(cyt, return_counts = True)[1]/(N))**2))
    Nall_null.append(NAll(nuc))
    Nall_infnull.append(len(np.unique(nuc[:,1])))
    for k in range(nMsats):
        List_AllFreqs_null.append(np.unique(nuc[:,k+2], return_counts = True)[1]/(N*2))
    Fis_null.append(Fis(N, nuc))
    Gen_null.append(Unique_genotype(N, nuc))
    Dist_mf_null.append(Dist_mf(nuc, False))
Dist_mf_null = np.array(Dist_mf_null)


# Stats Ne
print("-------------------------Ne-------------------------------")
Ne_trio = np.array(Ne_trio)
Ne_null = np.array(Ne_null)
print("thetaF Ne trio : ", np.mean(Ne_trio[:, 0]), "("+str(np.std(Ne_trio[:, 0]))+")")
print(mannwhitneyu(Ne_trio[:, 0], Ne_null[:, 0]))

print("thetaV Ne trio : ", np.mean(Ne_trio[:, 1]), "("+str(np.std(Ne_trio[:, 1]))+")")
print(mannwhitneyu(Ne_trio[:, 1], Ne_null[:, 1]))


print("He_trio :", np.mean(He_trio), "("+str(np.std(He_trio))+")")
H = np.mean(He_trio)
print("Ne trio estimation with inf allele locus", H/((1-H)*4*Muiloc))
print(mannwhitneyu(He_trio, He_null))


print("For cytoplamsic locus : ")
H = np.mean(Hec_trio)
print("trio mean H", H, "std : ", np.std(Hec_trio), "eff size", H/((1-H)*4*Muiloc))
H = np.mean(Hec_null)
print("null mean H", H, "std : ", np.std(Hec_null), "eff size", H/((1-H)*4*Muiloc))

# Stats Fis
print("-------------------------Fis------------------------------")
print("Fis trio : ", np.mean(Fis_trio), "Fis_null :",np.mean(Fis_null))
print(mannwhitneyu(Fis_trio, Fis_null))
print("median Fis !=0 : ", mannwhitneyu(Fis_trio, 0))
# Stats Nball
print("-------------------------Nall------------------------------")
print("Nall trio : ", np.mean(Nall_trio), np.std(Nall_trio), "Nall null :",np.mean(Nall_null), np.std(Nall_trio))
print(mannwhitneyu(Nall_trio, Nall_null))

print("inf all locus")
print("Nall trio : ", np.mean(Nall_inftrio), np.std(Nall_inftrio), "Nall null :",np.mean(Nall_infnull), np.std(Nall_infnull))
print("difference between mean at inf loc", np.mean(Nall_inftrio) -  np.mean(Nall_infnull))
print(mannwhitneyu(Nall_inftrio, Nall_infnull))


### Figs ####

#thetaF
plt.figure("thetaF")
bins = np.histogram_bin_edges(np.concat((Ne_trio[:,0], Ne_null[:,0])), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_null = np.histogram(Ne_null[:,0], bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

counts_trio = np.histogram(Ne_trio[:,0], bins=bins)[0]
counts_trio = counts_trio/np.sum(counts_trio)

plt.bar(x, counts_trio, w, alpha=0.7, label="Trioecy model ("+ str(non_fix)+" replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Ne estimated from microsat homozigoty")
plt.legend()
plt.show()

#thetaV
plt.figure("thetaV")
bins = np.histogram_bin_edges(np.concat((Ne_trio[:,1], Ne_null[:,1])), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_null = np.histogram(Ne_null[:,1], bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

counts_trio = np.histogram(Ne_trio[:,1], bins=bins)[0]
counts_trio = counts_trio/np.sum(counts_trio)

plt.bar(x, counts_trio, w, alpha=0.7, label="Trioecy model ("+ str(non_fix)+" replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Ne estimated from microsat variance in size")
plt.legend()
plt.show()


# Nb all
plt.figure("NbAll")
bins = np.histogram_bin_edges(np.concat((Nall_null, Nall_trio)), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_trio = np.histogram(Nall_trio, bins=bins)[0]
counts_trio = counts_trio/np.sum(counts_trio)

counts_null = np.histogram(Nall_null, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

plt.bar(x, counts_trio, w, alpha=0.7, label="Trioecy model ("+ str(non_fix)+" replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Mean number of allele")
plt.legend()
plt.show()

# plt.figure("NbAll per sex")
# Nball_sex = np.array(Nball_sex)
# Nball_sex = Nball_sex[np.all(np.logical_not(np.isnan(Nball_sex)), axis = 1)]
# plt.boxplot(np.array(Nball_sex), tick_labels = ["Herm", "Males","Males+CMS", "Females" ])
# plt.title("Mean allele number per loci whith size correction")
# plt.show()



### SFS

List_AllFreqs_trio = np.concat(List_AllFreqs_trio)
List_AllFreqs_null = np.concat(List_AllFreqs_null)

list_global = np.concat((List_AllFreqs_trio, List_AllFreqs_null))
plt.rcParams.update({'font.size': 18})

plt.figure("SFS")
plt.subplots_adjust(left=0.05, right=0.95, top=1, bottom=0.15)
plt.axes().set_ylim(0, 0.58)
bins = np.histogram_bin_edges(list_global, bins=np.linspace(0,1,11))
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]


counts_trio = np.histogram(List_AllFreqs_trio, bins=bins)[0]
counts_trio = counts_trio/np.sum(counts_trio)

counts_null = np.histogram(List_AllFreqs_null, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

bar_width = w * 0.4
plt.bar(x - bar_width/2, counts_trio, bar_width,
        label="Trioécie ("+ str(non_fix)+" réplicats)", color="red")

plt.bar(x + bar_width/2, counts_null, bar_width,
        label="Sans biais (100 réplicats)", color="grey")

labels = [f"{bins[i]:.1f}-{bins[i+1]:.1f}" for i in range(len(bins)-1)]

plt.xticks(x, labels, rotation=45)

plt.xlabel("Fréquences alléliques")
plt.ylabel("Fréquences")
plt.legend()
plt.show()

### Fis
plt.figure("Fis")

bins = np.histogram_bin_edges(np.concat((Fis_null, Fis_trio)), bins="auto")
x = (bins[0:-1:1] + bins[1::1])/2
w = bins[1] - bins[0]

counts_trio = np.histogram(Fis_trio, bins=bins)[0]
counts_trio = counts_trio/np.sum(counts_trio)

counts_null = np.histogram(Fis_null, bins=bins)[0]
counts_null = counts_null/np.sum(counts_null)

plt.bar(x, counts_trio, w, alpha=0.7, label="Trioecy model ("+ str(non_fix)+" replicates)", color ="red")
plt.bar(x, counts_null, w,  alpha=0.7, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.title("Fis")
plt.legend()
plt.show()

plt.figure("Fis per sex")
plt.subplots_adjust(left=0.05, right=0.95, top=1, bottom=0.05)
Fis_sex = np.array(Fis_sex)
Fis_sex = Fis_sex[np.all(np.logical_not(np.isnan(Fis_sex)), axis = 1)]
plt.boxplot(np.array(Fis_sex), tick_labels = ["Hermaphrodites", "Mâles","Mâles+SMC", "Femelles" ])
plt.hlines(y=0, xmin=0, xmax=5, linestyle="--")
plt.ylabel("$F_{IS}$")

plt.show()

##### Fst between sexes classes ####
Dist_mat = np.array(Dist_mat)
mean_mat = np.mean(Dist_mat, axis=0)

Dist_Distrib = np.array(Dist_Distrib)
Dist_mean = np.mean(Dist_Distrib, axis=0)
mat_pvalue = np.zeros((4,4), float)
for i in range(1000):
    mat_pvalue += np.int8(Dist_mean[i,:,:] > mean_mat)
#### Genotypes ####

# Nb & prop unique & eveness
Gen_trio = np.array(Gen_trio)
Gen_null = np.array(Gen_null)

plt.figure("Nb genotypes")
plt.title("nb of unique genotypes")
counts, bins = np.histogram(np.concat((Gen_null[:,0], Gen_trio[:,0])))
plt.hist(Gen_trio[:,0], label="Trioecy model ("+ str(non_fix)+" replicates)", color ="red")
plt.hist(Gen_null[:,0], alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.legend()
plt.show()

plt.figure("proportion unique genotypes")
plt.title("proportion of unique genotypes")
counts, bins = np.histogram(np.concat((Gen_null[:,1], Gen_trio[:,1])))
plt.hist(Gen_trio[:,1],  bins = bins, label="Trioecy model ("+ str(non_fix)+" replicates)", color ="red")
plt.hist(Gen_null[:,1], bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.legend()
plt.show()

plt.figure("Clonal eveness")
plt.title("clonal eveness")
counts, bins = np.histogram(np.concat((Gen_null[:,2], Gen_trio[:,2])))
plt.hist(Gen_trio[:,2],  bins = bins, label="Trioecy model ("+ str(non_fix)+" replicates)", color ="red")
plt.hist(Gen_null[:,2], bins = bins, alpha=0.5, label="Null model with equivalent Ne (100 replicates)", color="grey")
plt.legend()
plt.show()

plt.figure("Nb Genotype function of mean number of alleles")
plt.plot(Nall_trio/np.max(Nall_trio), Gen_trio[:,1], "o", color ="red")
plt.plot(Nall_null/np.max(Nall_null), Gen_null[:,1], "o", color ="grey")
plt.xlabel("Norm Mean allele number per loci")
plt.ylabel("N Genoytypes")
plt.show()

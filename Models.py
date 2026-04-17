import numpy as np
from BasicFunctions import *

##### Init Functions ######

def init_m1m2(N, N2, nMsats):
    """Intialization function for the first and second two separate sexes initial sex ratio at 1:1"""

    #initializing microsatellites each individual share the same haplotype at start
    #each individual is two line one for each chromosome first mother herited haplotype
    #the first column contain the information for the sex determination locus 0 : female recessive
    # and the second contain an neutral locus (infinite allele locus)
    nuc = np.zeros((N2, nMsats+2),dtype=int)
    #infinite allele model one neutral locus
    cyt = np.zeros((N),dtype=int)
    repeats = np.random.poisson(30, nMsats) + 5
    for i in range(N2):
        nuc[i,2:] = repeats
        if i<N and i%2==1:
            nuc[i*2+1,0] = 1
    return nuc, cyt

def OvulesPollenM3(N, nuc, cyt, HO, HP, MCMSP):
    """returns array containing the contribution of each individual to ovules and pollen"""
    Ovules, Pollen = np.zeros(N), np.zeros(N)
    all_males = nuc[1::2, 0]==1
    cms  = cyt[:, 0]==1

    #sex true false array
    males_cms = np.logical_and(all_males, cms)
    hermaphrodites = np.logical_and(np.logical_not(all_males), np.logical_not(cms))
    females = np.logical_and(np.logical_not(all_males), cms)
    males_nocms = np.logical_and(all_males, np.logical_not(cms))

    # filling the arrays
    Ovules[females] = 1
    Ovules[hermaphrodites] = HO
    Pollen[males_cms] = MCMSP
    Pollen[hermaphrodites] = HP
    Pollen[males_nocms] = 1
    return Ovules, Pollen

def init_CMS(N, N2, nMsats, HO, HP, MCMSP):
    """Intialization function for the trioecy model initialize each genotype at 1/4
    HO, HP, MCMS represent the relative weight of an individual to pollen/Ovule contribution"""

    #initializing microsatellites each individual share the same haplotype at start
    #each individual is two line one for each chromosome first mother herited haplotype
    #the first column contain the information for the sex determination locus 0 : female recessive
    # and the second contain an neutral locus (infinite allele locus)
    nuc = np.zeros((N2, nMsats+2),dtype=int)
    #infinite allele model one neutral locus and first locus for the CMS 1 (WT 0)
    cyt = np.zeros((N, 2),dtype=int)
    repeats = np.random.poisson(30, nMsats) + 5
    for i in range(N):
        nuc[i*2,2:] = repeats
        nuc[i*2+1,2:] = repeats
        #initializing each phenotype 4 * 1/4
        if i>=N//4 and i<N//2:
            # male plants
            nuc[i*2+1,0] = 1
        elif i>=N//2 and i<(N//4)*3:
            # female plants
            cyt[i, 0] = 1
        elif i>=(N//4)*3:
            # male + CMS plants
            nuc[i*2+1,0] = 1
            cyt[i, 0] = 1
    Ovules, Pollen = OvulesPollenM3(N, nuc, cyt, HO, HP, MCMSP)
    return nuc, cyt, Pollen, Ovules

######### Life cycle function ####################

def sex_m1(param, N, N2, Nsurv, nMsats, Musat, Muiloc, mut_ID, nuc, cyt, tmp_nuc, tmp_cyt, Ovules, Pollen):
    """sexual reproduction and mutations with no biais. param"""

    for i in range(Nsurv, N):
        #mate choice
        p1_id=np.random.choice(N, p = Ovules/np.sum(Ovules))*2
        p2_id=np.random.choice(N, p = Pollen/np.sum(Pollen))*2

        #recombination
        tmp_nuc[i*2,:] = recombi_ind(nuc[(p1_id,p1_id+1),], nMsats)
        tmp_nuc[i*2+1,:] = recombi_ind(nuc[(p2_id,p2_id+1),], nMsats)
        rnum = np.random.rand(2)

        #cytotype transmission
        tmp_cyt[i] = cyt[p1_id//2]


    #mutation of all offspring
    Nsurv2 = 2*Nsurv
    mutation_msat(N2 - Nsurv2, nMsats, Musat, tmp_nuc[Nsurv2:,2:])
    mutations_iloc(N2 - Nsurv2, Muiloc, tmp_nuc[Nsurv2:,1], mut_ID, 0)
    #mutate all the cytotypes
    mutations_iloc(N - Nsurv, Muiloc, tmp_cyt[Nsurv:], mut_ID, 1)
    return 1 - np.copy(tmp_nuc[1::2, 0]), np.copy(tmp_nuc[1::2, 0])


def sex_m2(param, N, N2, Nsurv, nMsats, Musat, Muiloc, mut_ID, nuc, cyt, tmp_nuc, tmp_cyt, Ovules, Pollen):
    """sexual reproduction and mutations with differential mortality. param = [Sm]"""

    Sm = param[0]
    i=Nsurv

    while i<N:
        #mate choice
        p1_id = np.random.choice(N, p = Ovules/np.sum(Ovules))*2
        p2_id = np.random.choice(N, p = Pollen/np.sum(Pollen))*2

        #recombination
        male_gam = recombi_ind(nuc[(p2_id,p2_id+1),], nMsats)
        if male_gam[0]!=1 or np.random.rand()<Sm:
            tmp_nuc[i*2,:] = recombi_ind(nuc[(p1_id,p1_id+1),], nMsats)
            tmp_nuc[i*2+1,:] = male_gam
            tmp_cyt[i] = cyt[p1_id//2]
            i+=1
    #mutation of all offspring
    Nsurv2 = 2*Nsurv
    mutation_msat(N2 - Nsurv2, nMsats, Musat, tmp_nuc[Nsurv2:,2:])
    mutations_iloc(N2 - Nsurv2, Muiloc, tmp_nuc[Nsurv2:,1], mut_ID, 0)
    #mutate all the cytotypes
    mutations_iloc(N - Nsurv, Muiloc, tmp_cyt[Nsurv:], mut_ID, 1)
    return 1 - np.copy(tmp_nuc[1::2, 0]), np.copy(tmp_nuc[1::2, 0])

def sex_CMS(param, N, N2, Nsurv, nMsats, Musat, Muiloc, mut_ID, nuc, cyt, tmp_nuc, tmp_cyt,  Ovules, Pollen):
    """Sexual reproduction and mutations for trioecy + CMS. param = [HO, HP, MCMSP, s, d]"""

    HO, HP, MCMSP, s, d = param

    i=Nsurv
    while i<N:
        aborted = False
        #Mother choice
        p1_id = np.random.choice(N, p = Ovules/np.sum(Ovules))
        # hermphrodite has a chance of selfing
        if cyt[p1_id,0]==0 and np.random.rand()<=s:
            if np.random.rand()<=d:
                aborted = True
            else:
                p2_id = p1_id
        else: #father choice
            p2_id = np.random.choice(N, p = Pollen/np.sum(Pollen))

        if not(aborted):
            #recombination
            tmp_nuc[i*2,:] = recombi_ind(nuc[(p1_id*2,p1_id*2+1),], nMsats)
            tmp_nuc[i*2+1,:] = recombi_ind(nuc[(p2_id*2,p2_id*2+1),], nMsats)

            #cytotype transmission
            tmp_cyt[i,:] = cyt[p1_id, :]
        i+=1

    #mutation of all offspring
    Nsurv2 = 2*Nsurv
    mutation_msat(N2 - Nsurv2, nMsats, Musat, tmp_nuc[Nsurv2:,2:])
    mutations_iloc(N2 - Nsurv2, Muiloc, tmp_nuc[Nsurv2:,1], mut_ID, 0)
    #mutate all the cytotypes
    mutations_iloc(N - Nsurv, Muiloc, tmp_cyt[Nsurv:, 1], mut_ID, 1)
    return OvulesPollenM3(N, tmp_nuc, tmp_cyt, HO, HP, MCMSP)



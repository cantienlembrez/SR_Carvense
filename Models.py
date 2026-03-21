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
    Pollen = [[-1],[0]] # first line containt id; second cumulative weight to pollen/Ovules contribution of individuals
    Ovules = [[-1],[0]]
    O_ID, P_ID = 0, 0
    for i in range(N):
        nuc[i*2,2:] = repeats
        nuc[i*2+1,2:] = repeats
        #initializing each phenotype 4 * 1/4
        if i<N//4:
            # hermaphrodite plants
            Ovules[0].append(i)
            Ovules[1].append(Ovules[1][O_ID] + HO)
            Pollen[0].append(i)
            Pollen[1].append(Pollen[1][P_ID] + HP)
            O_ID+=1
            P_ID+=1
        elif i>=N//4 and i<N//2:
            nuc[i*2+1,0] = 1
            # male plants
            Pollen[0].append(i)
            Pollen[1].append(Pollen[1][P_ID] + 1)
            P_ID+=1
        elif i>=N//2 and i<(N//4)*3:
            # female plants
            cyt[i, 0] = 1
            Ovules[0].append(i)
            Ovules[1].append(Ovules[1][O_ID] + 1)
            O_ID+=1
        else:
            # male + CMS plants
            nuc[i*2+1,0] = 1
            cyt[i, 0] = 1
            Pollen[0].append(i)
            Pollen[1].append(Pollen[1][P_ID] + MCMSP)
            P_ID+=1
    return nuc, cyt, Pollen, Ovules

######### Life cycle function ####################

def gen_incr(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt):
    """make one life cycle stop after offspring birth"""

    for i in range(N):
        #mate choice
        p1_id=np.random.randint(N)*2
        p2_id=np.random.randint(N)*2
        while (nuc[p2_id+1,0]==nuc[p1_id+1,0]):
            p2_id=np.random.randint(N)*2
        #passing mother as parent 1
        if nuc[p1_id+1,0]==1:
            p1_id, p2_id = p2_id, p1_id

        #recombination
        tmp_nuc[i*2,:] = recombi_ind(nuc[(p1_id,p1_id+1),], nMsats)
        tmp_nuc[i*2+1,:] = recombi_ind(nuc[(p2_id,p2_id+1),], nMsats)

        #cytotype transmission
        tmp_cyt[i] = cyt[p1_id//2]

    #mutation of all offspring
    mutation_msat(N2, nMsats, Musat, tmp_nuc[:,2:])
    mutations_iloc(N2, Muiloc, tmp_nuc[:,1])
    #mutate all the cytotypes
    mutations_iloc(N, Muiloc, tmp_cyt)
    return


def gen_incr_m2(N, N2, nMsats, Musat, Muiloc, nuc, cyt, tmp_nuc, tmp_cyt, Sm):
    """make one life cycle stop after offspring birth"""

    i=0
    while i<N:
        #mate choice
        p1_id=np.random.randint(N)*2
        p2_id=np.random.randint(N)*2
        while (nuc[p2_id+1,0]==nuc[p1_id+1,0]):
            p2_id=np.random.randint(N)*2
        #passing mother as parent 1
        if nuc[p1_id+1,0]==1:
            p1_id, p2_id = p2_id, p1_id

        #recombination
        Pollen = recombi_ind(nuc[(p2_id,p2_id+1),], nMsats)
        if Pollen[0]!=1 or np.random.rand()<Sm:
            tmp_nuc[i*2,:] = recombi_ind(nuc[(p1_id,p1_id+1),], nMsats)
            tmp_nuc[i*2+1,:] = Pollen
            tmp_cyt[i] = cyt[p1_id//2]
            i+=1
    #mutation of all offspring
    mutation_msat(N2, nMsats, Musat, tmp_nuc[:,2:])
    mutations_iloc(N2, Muiloc, tmp_nuc[:,1])
    #mutate all the cytotypes
    mutations_iloc(N, Muiloc, tmp_cyt)
    return

def choiceFromTab(tab, tablen, cumul):
    rnum = np.random.uniform(low = 0, high = cumul)
    Imin, Imax = 0, tablen
    interval = Imax - Imin
    while interval!=1:
        Iint = Imin + interval//2
        if rnum>tab[1][Iint]:
            Imin = Iint
        else:
            Imax = Iint
        interval = Imax - Imin
    return tab[0][Imax]

def gen_incr_CMS(N, N2, nMsats, Musat, Muiloc, HO, HP, MCMSP, s, d, nuc, cyt, Ovules, Pollen, tmp_nuc, tmp_cyt):
    """make one life cycle stop after offspring birth"""

    tmp_Ovules = [[-1],[0]]
    tmp_Pollen = [[-1],[0]]
    O_ID = 0
    P_ID = 0

    OvLen = len(Ovules[0])-1
    PoLen = len(Pollen[0])-1
    OvCum = Ovules[1][OvLen]
    PoCum = Pollen[1][PoLen]

    i=0
    while i<N:
        aborted = False
        #Ovule choice
        p1_id = choiceFromTab(Ovules, OvLen, OvCum)
        # hermphrodite has a chance of selfing
        if cyt[p1_id,0]==0 and np.random.rand()<=s:
            if np.random.rand()<=d:
                aborted=True
            else:
                p2_id = p1_id
        else:
            p2_id = choiceFromTab(Pollen, PoLen, PoCum)

        if not(aborted):
            #recombination
            tmp_nuc[i*2,:] = recombi_ind(nuc[(p1_id*2,p1_id*2+1),], nMsats)
            tmp_nuc[i*2+1,:] = recombi_ind(nuc[(p2_id*2,p2_id*2+1),], nMsats)

            #cytotype transmission
            tmp_cyt[i,:] = cyt[p1_id, :]

            # Attributing weight
            if tmp_nuc[i*2+1,0]==0 and cyt[p1_id,0]==0:
                #hermaphrodite
                tmp_Ovules[0].append(i)
                tmp_Ovules[1].append(tmp_Ovules[1][O_ID] + HO)
                tmp_Pollen[0].append(i)
                tmp_Pollen[1].append(tmp_Pollen[1][P_ID] + HP)
                O_ID+=1
                P_ID+=1
            elif  tmp_nuc[i*2+1,0]==1 and cyt[p1_id,0]==0:
                #male
                tmp_Pollen[0].append(i)
                tmp_Pollen[1].append(tmp_Pollen[1][P_ID] + 1)
                P_ID+=1
            elif tmp_nuc[i*2+1,0]==0 and cyt[p1_id,0]==1:
                #female
                tmp_Ovules[0].append(i)
                tmp_Ovules[1].append(tmp_Ovules[1][O_ID] + 1)
                O_ID+=1
            else:
                #male + CMS
                tmp_Pollen[0].append(i)
                tmp_Pollen[1].append(tmp_Pollen[1][P_ID] + MCMSP)
                P_ID+=1
            i+=1

    #mutation of all offspring
    mutation_msat(N2, nMsats, Musat, tmp_nuc[:,2:])
    mutations_iloc(N2, Muiloc, tmp_nuc[:,1])
    #mutate all the cytotypes
    mutations_iloc(N, Muiloc, tmp_cyt[:, 1])
    return tmp_Ovules, tmp_Pollen

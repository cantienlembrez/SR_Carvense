import numpy as np
from BasicFunctions import *

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
    mutations_iloc(N, Muiloc, tmp_cyt[:, 1])
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
    mutations_iloc(N, Muiloc, tmp_cyt[:, 1])
    return


def choiceFromTab(tab, tablen, cumul):
    rnum = np.random.uniform(low = 0, high = cumul)
    i=0
    while tab[1][i]<rnum:
        i+=1;
        #optimizable using minimal value and difference;
    return tab[0][i]


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

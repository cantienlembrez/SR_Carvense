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
    tmp_cyt = mutations_iloc(N, Muiloc, tmp_cyt)
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
        gam1 = recombi_ind(nuc[p2_id:(p2_id+1),:], nMsats)
        if gam1[0]!=1 or np.random.rand()<Sm:
            gam0 = recombi_ind(nuc[p1_id:(p1_id+1),:], nMsats)
            tmp_cyt[i] = cyt[p1_id//2]
            i+=1
    #mutation of all offspring
    mutation_msat(N2, nMsats, Musat, tmp_nuc[:,2:])
    mutations_iloc(N2, Muiloc, tmp_nuc[:,1])
    #mutate all the cytotypes
    tmp_cyt = mutations_iloc(N, Muiloc, tmp_cyt)
    return

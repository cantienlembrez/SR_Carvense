import numpy as np
import matplotlib.pyplot as plt

def FreqsAtEq(g, a, s, d, em):
    #See Nguyen and Pannell 2024 supplementary information
    if a < 2*(1-s*d)/(1-s):
        return None, None, None, None
    Px = (s*(1-d))/(g+s-1)
    Py = 1 - Px
    fXXn = (a*(2*Px-1)*(1-em)*g*Px) / (2*g*Px - a*(2*Px-1)*(1-s-(1-em)*(1-s*d)))
    fXYn = (fXXn * (1-s)*Py) / (g*Px)
    fXXc = Px - (fXXn*(1-s*d)) / g
    fXYc = (fXXc*Py) / Px
    if any([fXXn<=0, fXYn+fXYc<=0, fXXc<=0]):
        return None, None, None, None
    return fXXn, fXYn, fXXc, fXYc


d = 0
EM = np.linspace(0.1, 0.9, 3)
S = np.linspace(0.3, 0.5, 3)
G = np.linspace(0.7, 2, 100)
A = np.linspace(1.5, 14, 100)

plt.figure()

f, axarr = plt.subplots(len(S),len(EM))

for l in range(len(EM)):
    em = EM[l]
    for c in range(len(S)):
        mat = np.zeros((len(G),len(A)), dtype=float)
        s = S[c]
        for ig in range(len(G)):
            g = G[ig]
            for ia in range(len(A)):
                a = A[ia]
                mat[len(G)-ig-1, ia] = FreqsAtEq(g, a, s, d, em)[0]
        m = axarr[c,l].imshow(mat, extent=[A[0], A[len(A)-1], G[0], G[len(G)-1]], aspect="auto", vmin=0, vmax=1)
        axarr[c,l].set_title("s : "+str(round(s,1))+" em : "+str(round(em, 1)))
        axarr[c,l].set_ylim(G[0], G[len(G)-1])
f.colorbar(m)
plt.show()

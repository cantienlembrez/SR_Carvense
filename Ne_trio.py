from sympy import *

init_printing()

# Determining effective pop size of Nguyen and Pannell 2025 model
# Using Charlesworth and Charlesworth 2010 method

H, M, C, F = symbols("N_h N_m N_c N_f")
alpha, g, em, s, d = symbols("a g e_m s d")
Fh = symbols("Fh")
# ClassContribution to each classes for beta coefficients
#prs = p(r <- s)

# contrib to M
pmf = Integer(0)
pmh = Rational(1, 2)
pmm = Rational(1, 2)*((alpha*M)/(alpha*M + (1-em)*alpha*C))
pmc = Rational(1, 2)*((alpha*(1-em)*C)/(alpha*M + (1-em)*alpha*C))

print(simplify(pmf+pmh+pmm+pmc))

# contrib to C
pcf = Rational(1, 2)
pch = Integer(0)
pcm = Rational(1, 2)*((alpha*M)/(alpha*M + (1-em)*alpha*C))
pcc = Rational(1, 2)*((alpha*(1-em)*C)/(alpha*M + (1-em)*alpha*C))

print(simplify(pcf+pch+pcm+pcc))

# contrib to F
pff = Rational(1, 2)
pfh = Rational(1, 2)*((2*H)/(alpha*M + (1-em)*alpha*C + 2*H))
pfm = Rational(1, 2)*((alpha*M) /(alpha*M + (1-em)*alpha*C + 2*H))
pfc = Rational(1, 2)*(((1-em)*alpha*C) /(alpha*M + (1-em)*alpha*C + 2*H))

print(simplify(pff+pfh+pfm+pfc))

#contrib to H
seff = s*(1-d)/(1-s*d)
phf = Integer(0)
phh = Rational(1, 2) + Rational(1, 2)*(seff + (1-seff)*((2*H)/(alpha*M + (1-em)*alpha*C + 2*H)))
phm = Rational(1, 2)*(1-seff)*((alpha*M )/(alpha*M + (1-em)*alpha*C + 2*H))
phc = Rational(1, 2)*(1-seff)*(((1-em)*alpha*C)/(alpha*M + (1-em)*alpha*C + 2*H))

print(simplify(phf+phh+phm+phc))

# tab of contribution for computing the beta coefficents
# lines : offspring and cols : parent
#   H   M   C   F
contrib_tab = [
  [phh,phm,phc,phf], #H
  [pmh,pmm,pmc,pmf], #M
  [pch,pcm,pcc,pcf], #C
  [pfh,pfm,pfc,pff]  #F
]

flux_matrix = Matrix(contrib_tab)

# alpha_i are deducted from flux matrix

a = flux_matrix.T.eigenvects()
v = a[1][2][0] #1 is always the dominant eig ?
a = v/sum(v)
# theta coeficients theta_rsu = 1/N_u for all (r,s)
# H M C F
theta = [
 1/H, 1/M, 1/C, 1/F
]

# gamma coefficient 1/2 for all u!=H for u=H 1/2+1/2*(Fis)
gamma = [
 Rational(1,2)+Rational(1,2)*Fh, Rational(1,2), Rational(1,2), Rational(1,2)
]
Idx = {
    "H":0,"M":1,"C":2,"F":3
}


Pc = Rational(0)
for R in ["H", "M", "C", "F"]:
    for S in ["H", "M", "C", "F"]:
        for U in ["H", "M", "C", "F"]:
            beta = contrib_tab[Idx[R]][Idx[U]] * contrib_tab[Idx[S]][Idx[U]]
            Pc += a[Idx[R]] * a[Idx[S]] * beta * theta[Idx[U]] * gamma[Idx[U]]



substitutions = {
 g:1.1, alpha:4.8, em:0.9, s:0.8, d:0.3,
 Fh:0.68661903, H:65.64364876, M: 7.24637681, C:350.24154589, F:576.86842853
 }


# substitutions = {
# g:1.1, alpha:4.8, em:0.9, s:0.8, d:0.3,
# Fh:6.58661903e-01, H:64.43010753, M:6.78494623655914, C:355.35483871, F:573.43010753
# }


# substitutions = {
# g:1.1, alpha:4.8, em:0.9, s:0.8, d:0.3,
# Fh:0.68570668, H:330.17, M:36.47, C:1751.18, F:2882.18
# }

Pc_num = Pc.subs(substitutions)
print("Ne", 1/(2*Pc_num))
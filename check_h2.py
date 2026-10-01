import math
import itertools

exp3 = [3.42525091, 0.62391373, 0.16885540]
c = [0.15432897, 0.53532814, 0.44463454]
def nrm(a): return (2*a/math.pi)**0.75
R = 1.4
RA, RB = -R/2, R/2

def ovlp_raw(i, j, CA, CB):
    p = exp3[i]+exp3[j]
    Kab = math.exp(-exp3[i]*exp3[j]/p*(CA-CB)**2)
    return Kab*(math.pi/p)**1.5

def E_full(i, j, k, l):
    cen = {0: RA, 1: RB}
    p = exp3[i]+exp3[j]; q = exp3[k]+exp3[l]
    P = (exp3[i]*cen[i]+exp3[j]*cen[j])/p
    Q = (exp3[k]*cen[k]+exp3[l]*cen[l])/q
    Kab = math.exp(-exp3[i]*exp3[j]/p*(cen[i]-cen[j])**2)
    Kcd = math.exp(-exp3[k]*exp3[l]/q*(cen[k]-cen[l])**2)
    T = p*q/(p+q)*(P-Q)**2
    F0 = 0.5*math.sqrt(math.pi/T)*math.erf(math.sqrt(T)) if T > 1e-12 else 1.0
    nn = nrm(exp3[i])*nrm(exp3[j])*nrm(exp3[k])*nrm(exp3[l])
    cc = c[i]*c[j]*c[k]*c[l]
    return 2*math.pi**2.5/(p*q*math.sqrt(p+q))*Kab*Kcd*F0*nn*cc

def one_electron(fn, C1, C2):
    # normalized contracted function at C1 vs C2: includes renorm factor f^2
    s_self = 0.0
    for i in range(3):
        for j in range(3):
            s_self += c[i]*c[j]*nrm(exp3[i])*nrm(exp3[j])*ovlp_raw(i, j, C1, C1)
    f2 = 1.0/s_self
    s = 0.0
    for i in range(3):
        for j in range(3):
            s += c[i]*c[j]*nrm(exp3[i])*nrm(exp3[j])*fn(i, j, C1, C2)
    return s*f2

def kinetic_raw(i, j, C1, C2):
    return 3*exp3[i]*exp3[j]/(exp3[i]+exp3[j])*ovlp_raw(i, j, C1, C2)

def nuclear_raw(i, j, C1, C2, CC):
    p = exp3[i]+exp3[j]
    P = (exp3[i]*C1+exp3[j]*C2)/p
    T = p*(P-CC)**2
    F0 = 0.5*math.sqrt(math.pi/T)*math.erf(math.sqrt(T)) if T > 1e-12 else 1.0
    return -2*math.pi/p*ovlp_raw(i, j, C1, C2)*F0

h11 = one_electron(lambda i, j, A, B: kinetic_raw(i, j, A, B) + nuclear_raw(i, j, A, B, A), RA, RA)
h12 = one_electron(lambda i, j, A, B: kinetic_raw(i, j, A, B) + nuclear_raw(i, j, A, B, A), RA, RB)
s12 = one_electron(lambda i, j, A, B: ovlp_raw(i, j, A, B), RA, RB)

def eri(i, j, k, l):
    cen = {0: RA, 1: RB}
    s_self = 0.0
    for a in range(3):
        for b in range(3):
            s_self += c[a]*c[b]*nrm(exp3[a])*nrm(exp3[b])*ovlp_raw(a, b, RA, RA)
    f4 = 1.0/s_self**2
    p = exp3[i]+exp3[j]; q = exp3[k]+exp3[l]
    P = (exp3[i]*cen[i]+exp3[j]*cen[j])/p
    Q = (exp3[k]*cen[k]+exp3[l]*cen[l])/q
    Kab = math.exp(-exp3[i]*exp3[j]/p*(cen[i]-cen[j])**2)
    Kcd = math.exp(-exp3[k]*exp3[l]/q*(cen[k]-cen[l])**2)
    T = p*q/(p+q)*(P-Q)**2
    F0 = 0.5*math.sqrt(math.pi/T)*math.erf(math.sqrt(T)) if T > 1e-12 else 1.0
    nn = nrm(exp3[i])*nrm(exp3[j])*nrm(exp3[k])*nrm(exp3[l])
    cc = c[i]*c[j]*c[k]*c[l]
    return 2*math.pi**2.5/(p*q*math.sqrt(p+q))*Kab*Kcd*F0*nn*cc*f4

print("h11 =", round(h11, 6), " h12 =", round(h12, 6), " s12 =", round(s12, 6))
print("(11|11) =", round(eri(0, 0, 0, 0), 5))
print("(12|12) =", round(eri(1, 1, 1, 1) if False else eri(0, 1, 0, 1), 5))
print("(11|22) =", round(eri(0, 0, 1, 1), 5))
print("(11|12) =", round(eri(0, 0, 0, 1), 5))

# SCF (symmetric MO is the exact solution)
w2 = 1.0/(2*(1+s12))
D = [[2*w2]*2 for _ in range(2)]
H = [[h11, h12], [h12, h11]]
def fock(i, j):
    g = 0.0
    for k in range(2):
        for l in range(2):
            g += D[k][l]*(eri(i, j, k, l) - 0.5*eri(i, k, j, l))
    return H[i][j] + g
E = 0.0
for mu in range(2):
    for nu in range(2):
        E += D[mu][nu]*(H[mu][nu] + fock(mu, nu))
E *= 0.5
print("E_elec =", round(E, 8))
print("E_total =", round(E + 1.0/R, 8), " (literature -1.11675931)")

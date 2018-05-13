import scipy.misc
import numpy
# eerste stap: grootheden


def density(radius):
    import math.pi as pi
    ro = 1/((2*pi)*radius*(radius+1)**3)
    return ro


def mass(radius):
    r = radius
    m = (r**2)/(1+r)**2
    return m


def BindPot(r):
    V = 1/(r+1)
    return V


def BindPotDer1(r):
    fprime = scipy.misc.derivative(BindPot, r)
    return fprime


def BindPotDer2(r):
    fprime = scipy.misc.derivative(BindPot, r, n=2)
    return fprime


# tweede stap: Routines om heen en weer te gaan tussen
# Energie E en draaimoment L, pericentrumafstand en
# apocentrumafstand


def aperi(E, L):

    # De eerste lijnen is het oplossen van de vergelijking uit 1.112
    p = [2*E, 2*E-2, L**2, L**2]
    Roots = numpy.roots(p)
    aperi = set()
    # De roots eruit filteren die complex zijn
    RRoots = set(Roots.real[abs(Roots.imag) < 1e-5])  # where I chose 1-e5 as a threshold
    # De roots eruit filteren die negatief zijn
    for root in RRoots:
        if root > -0.00001:
            aperi.add(round(root, 4))
    return aperi
# Volgende functies geven minimum en maximum oplossing van 1.112


def aphelium(E, L):
    aperium = list(aperi(E, L))
    return max(aperium)


def perihelium(E, L):
    aperium = list(aperi(E, L))
    return min(aperium)


def E(ap, peri):
    a = numpy.array([[2*(ap**3)+2*(ap**2), ap+1], [2*(peri**3)+2*(peri**2), peri + 1]])
    b = numpy.array([2*(ap**2), 2*(peri**2)])
    return(numpy.linalg.solve(a, b)[0])


def L(ap, peri):
    a = numpy.array([[2*(ap**3)+2*(ap**2), ap+1], [2*(peri**3)+2*(peri**2), peri + 1]])
    b = numpy.array([2*(ap**2), 2*(peri**2)])
    return(numpy.sqrt(numpy.linalg.solve(a, b)[1]))

# Radiele periode volledig bepaald door peri-en apoheleum
# want men kan via peri en apo de energie en draaimoment bepalen


def T_rad(peri, apo):
    import scipy.integrate as integrate
    E_ = E(apo, peri)
    L_ = L(apo, peri)
# integrate.quad neemt enkel functie objecten of methoden aan

    def functie(r):
        return 2 / (2*(E_ + BindPot(r)) - L_**2 / r**2)**0.5
    return integrate.quad(functie, peri, apo)[0]


# Stap 3: Rosettebanen integreren
# waiting for Triss to finish

# Stap 4: voor verschillende r'en de E en L berekenen
# Dan een fit maken zodat men voor willekeurige E de L weet

# een functie die de r geeft bij een gegeven massa, uitgedrukt in decimalen
# Totale massa is 1
# Dus bij 99% van de totale massa is n = 0.99
def r_mass(n):
    r_half = 1 + 2**0.5
    r_max = (((2*n)**0.5)*r_half) / (1 + (1 - (2*n)**0.5)*r_half)
    # bij M = 0.99 is dit ongeveer 198.5
    return r_max


def findL(E):
    r_max = r_mass(0.99)  # 99% van de totale massa wordt beschouwd
    E = []  # x-as
    L = []  # y-as
    r = 1
    for r in range(0, r_max):
        M = mass(r_max)
        # ingesloten massa, zou moeten gelijk zijn aan n zoals gedefinieerd
        # boven de functie r_mass
        orb_m = (M*r)**(0.5)
        e = BindPot(r) + M/(2*r)
        E.insert(e)[0]
        L.insert(orb_m)[0]
        # E moet in stijgende volgorde zijn voor de komende plotfunctie
        # E daalt bij hogere r, vandaar insert
    spl = scipy.interpolate.UnivariateSpline(E, L)
    return spl.__call__(E)  # returnt de waarde van de fit op positie E

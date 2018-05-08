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

# Volgende definities gelden enkel op E = V_eff

def aperi(E, L):
    #De eerste lijnen is het oplossen van de vergelijking uit 1.112
    p=[2*E, 2*E-2, L**2, L**2]
    Roots = numpy.roots(p)
    aperi=set()
    #De roots eruit filteren die complex zijn
    RRoots = set(Roots.real[abs(Roots.imag)<1e-5]) # where I chose 1-e5 as a threshold
    #De roots eruit filteren die negatief zijn
    for root in RRoots:
        if root > -0.00001:
            aperi.add(round(root,4))
    return aperi
#Volgende functies geven minimum en maximum oplossing van 1.112    
def aphelium(E,L):    
    aperium = list(aperi(E,L))
    return max(aperium)

def perihelium(E,L):
    aperium = list(aperi(E,L))
    return min(aperium)
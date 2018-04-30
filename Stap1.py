import scipy.misc
from math import sqrt
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

def EtoL(E, r):
    L = sqrt(BindPot(r) - 2*r**2*E)
    return L


def LtoE(L, r):
    E = (L**2)/(2*r**2) + BindPot(r)
    return E

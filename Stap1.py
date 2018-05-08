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
    """ De eerste lijnen is het oplossen van de vergelijking uit 1.112 """
<<<<<<< HEAD
    p = [2*E, 2*E-2, L**2, (L**2)]
    aperi = numpy.roots(p)
    i = 0
    """Hier filteren we alle oplossingen die positief zijn, en dus fysisch eruit"""
    while i < 3:
        if aperi[i] > 0:
            del aperi[i]
            i = 3
        i += 1
    """aperi is nu een tuple van alle fysische oplossingen"""
    return(aperi)
    """Vervolgens bepalen we welke aphelium en welke periheliumafstand is"""
    """Eventueel kan dit herschreven worden in 2 aparte functies, een voor peri een voor aphelium"""


def peri(E, L):
    aper = aperi(E, L)
    if aper[0] > aper[1]:
        return aper[1]
    else:
        return aper[0]

def ap(E, L):
    aper = aperi(E, L)
    if aper[0] > aper[1]:
        return aper[0]
    else:
        return aper[1]
=======
    p=[2*E, 2*E-2, L**2, (L**2)]
    aperi = numpy.roots(p)
    i = 0
    aperi2=[]
    """Hier filteren we alle oplossingen die positief zijn, en dus fysisch eruit"""
    while i<3:
        if aperi[i]>0:
            del aperi[i]
            i = 3
        i+=1
    return(aperi)
    """Vervolgens bepalen we welke aphelium en welke periheliumafstand is"""
    """Eventueel kan dit herschreven worden in 2 aparte functies, een voor peri een voor aphelium"""
    if aperi[0]>aperi[1]:
        aphelium = aperi[0]
        perihelium = aperi[1]
    elif aperi[0]==aperi[1]:
        aphelium = aperi[1]
        perihelium = aperi[0] 
    else:
        aphelium = aperi2[1]
        perihelium = aperi2[0]
    return(perihelium, aphelium)
>>>>>>> 0d101ae349d3c751746bef29d75a9c03ac8e079c

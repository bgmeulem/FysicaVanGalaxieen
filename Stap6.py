# Stap 6: berkenenen welke massa we aan sterren moeten geven om het
# Hernquistmodel te reproduceren

from Stap5 import rad_distr_e, mass_increase
import numpy as np
from scipy.optimize import leastsq
import time


def BijdrageL(E, mass_frac=0.9, rinterval=100):
    l_distr = rad_distr_e(mass_frac, E, rinterval)
    # een lijst van alle radiële distributies (elk afhankelijk van een L)
    # die bij 1 E horen
    # Deze functie sommeert alle elementen op dezelfde index en steekt ze
    # in een nieuwe lijst
    bijdrage = [float(sum(elements))/float(len(l_distr)) for elements in zip(*l_distr)]
    return bijdrage
    # Dit is f_i(E_k)
    # oftewel de bijdrage van 1 E (en alle L) op een stukje r_i


def m_k(massfrac=0.9, rinterval=100):
    MatrixE = []
    for E in np.linspace(0.0001, 0.999, rinterval):
        Edistr = (BijdrageL(E, massfrac, rinterval))
        MatrixE.append(list(Edistr))

    # MatrixE is nu van de vorm [[Alle bijdrages van E1 voor 0-rmax],[E2],...]
    v = np.array([list(x) for x in MatrixE])
    # v is nu een array van lijsten van lijsten van de vorm [[Alle bijdrages van E1 voor 0-rmax],[E2],...]
    M = (mass_increase(massfrac, rinterval))
    # Dit is de theoretisch bepaalde Massatoenames tussen de verschillende intervallen
    M = np.asarray(M)
    # numpy.dot werkt liever met arrays
    mk_init = np.array(list(1.0/rinterval for i in range(rinterval)))
    # Een uniforme verdeling als initiele schatting voor mk
    
    # hier gebeurt de optimalisatie via leastsquare
    def f(mk, M, v):
        return M - np.dot(mk, v)
    m_optimal = leastsq(f, mk_init, args=(M, v))[0]
    return(m_optimal)

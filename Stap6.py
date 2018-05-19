# Stap 6: berkenenen welke massa we aan sterren moeten geven om het
# Hernquistmodel te reproduceren

from Stap5 import rad_distr_e, mass_increase
import numpy as np
from Stap4 import r_mass
from scipy.optimize import leastsq


def BijdrageL(E, mass_frac, i):
    l_distr = rad_distr_e(mass_frac, E, i)
    # een lijst van alle radiële distributies (elk afhankelijk van een L)
    # die bij 1 E horen
    # Deze functie sommeert alle elementen op dezelfde index en steekt ze
    # in een nieuwe lijst
    bijdrage = [float(sum(elements))/float(len(l_distr)) for elements in zip(*l_distr)]
    return bijdrage
    # Dit is f_i(E_k)
    # oftewel de bijdrage van 1 E (en alle L) op een stukje r_i


def m_k(massfrac, stapjes):
    stapjes = 20
    massfrac = 0.8
    MatrixE = []
    for E in np.linspace(0.4, 0.9, stapjes):
        Edistr = (BijdrageL(E, massfrac, stapjes))
        MatrixE.append(list(Edistr))
    v = [list(x) for x in zip(*MatrixE)]
    # MatrixE is nu van de vorm [[Alle bijdrages van E1 voor 0-rmax],[E2],...]
    # MatrixE is van de vorm [[Alle E-bijdrages voor r0], [Alle voor r1], ...]

    M = np.asarray(mass_increase(0.99))
    # uniforme verdelen als gok om te beginnen
    mk_init = np.asarray(list(1.0/stapjes for i in range(stapjes)))

    def f(mk, v, M):
        return M - np.dot(v, mk)

    m_optimal = leastsq(f, mk_init, args=(M, v))
    print(m_optimal)
    return(m_optimal)
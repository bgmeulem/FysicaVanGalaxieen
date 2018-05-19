# Stap 6: berkenenen welke massa we aan sterren moeten geven om het
# Hernquistmodel te reproduceren

from Stap5 import rad_distr_e, mass_increase
from numpy import linspace, dot
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
    massfrac = 0.99
    interval = linspace(0, r_mass(massfrac), stapjes)
    MatrixE = []
    for E in linspace(0.4, 0.9, stapjes):
        Edistr = BijdrageL(E, 0.99, stapjes)
        for element in Edistr:
            E_ri = []
            E_ri.append(element)
            # E_ri is nu de E-bijdrage van alle E's voor het interval ri
        MatrixE.append(E_ri)
    # MatrixE is van de vorm [[Alle E-bijdrages op r1], [Alle op r2], ...]

    def f(mk, M, v):
        return M - dot(v, mk)
    M = mass_increase(0.99)
    v = MatrixE
    mk_init = list(1.0/len(interval) for i in range(interval - 1))
    # uniforme verdelen als gok om te beginnen
    m_optimal = leastsq(f, mk_init, args=(M, v))
    print(m_optimal)
    return(m_optimal)

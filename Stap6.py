# Stap 6: berkenenen welke massa we aan sterren moeten geven om het
# Hernquistmodel te reproduceren

from Stap5 import rad_distr_e, mass_increase
import numpy as np
from scipy.optimize import leastsq
import time

begin = time.time()

def BijdrageL(E, mass_frac, rinterval):
    l_distr = rad_distr_e(mass_frac, E, rinterval)
    # een lijst van alle radiële distributies (elk afhankelijk van een L)
    # die bij 1 E horen
    # Deze functie sommeert alle elementen op dezelfde index en steekt ze
    # in een nieuwe lijst
    bijdrage = [float(sum(elements))/float(len(l_distr)) for elements in zip(*l_distr)]
    return bijdrage
    # Dit is f_i(E_k)
    # oftewel de bijdrage van 1 E (en alle L) op een stukje r_i


def m_k(massfrac, rinterval):
    MatrixE = []
    for E in np.linspace(0.0001, 0.999, rinterval):
        Edistr = (BijdrageL(E, massfrac, rinterval))
        MatrixE.append(list(Edistr))

    # MatrixE is nu van de vorm [[Alle bijdrages van E1 voor 0-rmax],[E2],...]
    v = np.array([list(x) for x in zip(*MatrixE)])
    # MatrixE is van de vorm [[Alle E-bijdrages voor r0], [Alle voor r1], ...]

    M = (mass_increase(0.99, rinterval))
    M.append(0)  # er komt nagenoeg 0 massa bij na r_max
    # dit is om de dimensies te doen kloppen voor de matrixbewerking
    # uniforme verdelen als gok om te beginnen
    M = np.asarray(M)
    # numpy.dot werkt liever met arrays
    mk_init = np.array(list(1.0/rinterval for i in range(rinterval)))
    # Een uniforme verdeling als initiele schatting voor mk?

    def f(mk, M, v):
        return M - np.dot(mk, v)
    m_optimal = leastsq(f, mk_init, args=(M, v))[0]
    return(m_optimal)


print(m_k(0.99, 1000))
end = time.time()
print(end-begin)

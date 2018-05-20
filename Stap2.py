# tweede stap: Routines om heen en weer te gaan tussen
# Energie E en draaimoment L, pericentrumafstand en
# apocentrumafstand
from numpy import sqrt
from numpy import roots
from Stap1 import BindPot
from numpy import pi


def aperi(E, L):
    aperi_list = []
    if L < 10**(-4):
        # L is praktisch 0: ster oscilleert of zit stil
        aperi_list.append(float(0.0))
        if E != 0:  # ster oscilleert
            aperi_list.append((1-E)/E)
        else:  # ster zit stil
            aperi_list.append(float(0.0))
    else:  # L != 0: cirkelbaan of ellips
        oplossingen = list(roots([2*E, 2*E - 2, L**2, L**2]))
        for element in oplossingen:
            # complexe en negatieve wortels filteren
            if isinstance(element, float) and element > 0:
                aperi_list.append(float(element))
    aperi_list.sort()
    return aperi_list

# Volgende functies geven minimum en maximum oplossing van 1.112


def aphelium(E, L):
    return aperi(E, L)[-1]


def perihelium(E, L):
    return aperi(E, L)[0]


def energie(ap, peri):
    if ap > 0:
        return((ap + peri + ap*peri)/((ap + peri)*(ap + 1)*(peri+1)))
    else:
        return 1


def draaimoment(ap, peri):
    if ap == 0:
        return 0
    return sqrt((2*(ap**2)*(peri**2))/((ap + peri)*(ap + 1)*(peri + 1)))

# Radiele periode volledig bepaald door peri-en apoheleum
# want men kan via peri en apo de energie en draaimoment bepalen


# in deze functie zit ook nog de periode voor een cirkelbaan, dit om
# de dingen wat te bundelen
def T_rad(apo, peri):

    import scipy.integrate as integrate
    E = energie(apo, peri)
    L = draaimoment(apo, peri)
    # integrate.quad neemt enkel functie objecten of methoden aan

    def functie(r):
        return 2*sqrt(1 / (2*(-E + BindPot(r)) - (L**2 / r**2)))

    if peri != apo:
        periode = abs(integrate.quad(functie, apo, peri)[0])
        return periode
    else:
        return (pi*2*(apo**2))/L
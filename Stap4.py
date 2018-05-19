# Stap 4: voor verschillende r'en de E en L berekenen
# Dan een fit maken zodat men voor willekeurige E de L weet

from scipy.interpolate import UnivariateSpline
from numpy import linspace
from Stap2 import draaimoment, energie

# een functie die de r geeft bij een gegeven massa, uitgedrukt in decimalen
# Totale massa is 1
# Dus bij 99,9% van de totale massa is n = 0.999


def r_mass(n):
    r_half = 1 + 2**0.5
    r_max = (((2*n)**0.5)*r_half) / (1 + (1 - (2*n)**0.5)*r_half)
    # bij M = 0.999 is dit ongeveer 1998.5
    return r_max

# 99.9% van de totale massa wordt beschouwd
r_max = r_mass(0.999)
E = []
L = []
for r in linspace(0, r_max, 100000):
    orb_mom = draaimoment(r, r)
    e = energie(r, r)
    E.insert(0, e)
    L.insert(0, orb_mom)
    # E moet in stijgende volgorde zijn voor de komende plotfunctie
    # E daalt bij hogere r, vandaar insert
    # returnt een lijst met als eerste element een lijst van alle E waarden
    # en als tweede de L-waarden
ListELcouples = [E,L]
# dit geldt enkel om het draaimoment te vinden van cirkelbanen
E_list = ListELcouples[0]
L_list = ListELcouples[1]
spl = UnivariateSpline(E_list, L_list, s=0)



def findL(E):
    # plt.plot(numpy.arange(0, 1, 0.0001), spl(numpy.arange(0, 1, 0.0001)))
    # plt.plot(E_list, L_list)
    # plot neemt kwadratisch af naar 0, zoals het hoort
    return spl.__call__(E)  # returnt de waarde van de fit op positie E
# Stap 5: interval 0:r_max opdelen in gelijke intervallen
# Massatoename berekenen
# Integraalruimte opdelen in grid (E_k, L_k)
from numpy import linspace, histogram
from Stap1 import mass
from Stap4 import r_mass, fitL, findL
from Stap3 import BaanInt
from Stap2 import aphelium, perihelium, T_rad
import matplotlib.pyplot as plt


def interval_r(r_max, n):
    return linspace(0, r_max, n)


def mass_increase(r1, r2):
    if r1:
        return mass(r2) - mass(r1)
    else:
        return mass(r1)


def rad_distr_e(r_max, e, i=100):
    # een radiele distributie voor 1 zekere e
    # het straal-interval wordt standaard verdeeld in 100 stukjes
    r_max = r_mass(0.4)
    interval = linspace(0, r_max, i)
    rad_distr_E = list()
    spl = fitL(e)
    for l in linspace(0, findL(e, spl), 20):
        L = findL(e, spl)
        # Draaimoment bij cirkelbaan is steeds de maximale voor een
        # bepaalde energie
        apo = aphelium(e, L)
        peri = perihelium(e, L)
        baan_rad = BaanInt(apo, peri, 1000)[1]
        baan_rad_half = baan_rad[:(len(baan_rad)//2)]
        # histogram normaliseren op de halve radiële periode (tot op orde -16)
        weights = len(baan_rad_half)*[1.0/len(baan_rad_half)]
        # histogram verdeelt de radiële distributie in bins, deze bins
        # worden bepaald door ons interval
        fractie = tuple(histogram(baan_rad_half, interval, weights=weights)[0])
        rad_distr_E.append(fractie)
        # plt.plot(tuple(histogram(baan_rad_half, interval, weights=weights)[1][:-1]), fractie)
        # bovenstaande plot kan getest worden voor 1 E-L koppel
        # run deze plot aub niet in iteratie met L'en en E's, dat gaat
        # poepeloeri lang duren
    # nu is rad_distr_E een lijst met als elk element de radiële distributie
    # (op hun beurt tuples) bij 1 L-waarde. Elk element van rad_distr_E
    # heeft dus 20 elementen (We verdelen in 20 L-waarden) en elk element
    # is een tuple van de radiële verdeling bij een E-L koppel
    return rad_distr_E


def rad_distr_tot(r_max, i):
    # i is het aantal delen dat we de r_max opdelen
    # een interval opgesteld van 0 tot r_max in 100 stukjes
    sterren_fractie = list()
    # nu gaan we rad_distr_e sommeren ove e om zo een verdeling van sterren
    # te krijgen
    for e in linspace(0, 0.999, 20):
        distr = list(rad_distr_e(r_max, e, i))
        sterren_fractie.append(distr)
    return sterren_fractie

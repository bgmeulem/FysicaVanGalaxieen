# Stap 5: interval 0:r_max opdelen in gelijke intervallen
# Massatoename berekenen
# Integraalruimte opdelen in grid (E_k, L_k)
from numpy import linspace, histogram
from Stap1 import mass
from Stap4 import r_mass, fitL, findL
from Stap3 import BaanInt
from Stap2 import aperi, aphelium, perihelium
# import matplotlib.pyplot as plt

spl = fitL()


def mass_increase(massfrac, i=100):
    r_max = r_mass(massfrac)
    interval = linspace(0, r_max, i)
    massincrease = []
    for n in range(len(interval) - 1):
        increase = mass(interval[n + 1]) - mass(interval[n])
        massincrease.append(increase)
    massincrease.append(0.0)  # er komt 0 massa bij tussen r_max en ??
    return massincrease
print(sum(mass_increase(0.9, 10)))

def rad_distr_e(mass_frac, e, i=100):
    # een radiele distributie voor 1 zekere e
    # het straal-interval wordt standaard verdeeld in 100 stukjes
    interval = linspace(0, r_mass(mass_frac), i+1)
    rad_distr_E = list()
    for l in linspace(findL(e, spl)/1000, findL(e, spl), 41):
        # Draaimoment bij cirkelbaan is steeds de maximale voor een
        # bepaalde energie
        if len(aperi(e, l)) != 2:  # geen fysische oplossing
            break
        apo = aphelium(e, l)
        peri = perihelium(e, l)
        baan_rad = BaanInt(apo, peri, 10000)[1]
        baan_rad_half = baan_rad[:(len(baan_rad)//2)]
        # histogram normaliseren op de halve radiele periode (tot op orde -16)
        weights = len(baan_rad_half)*[1.0/len(baan_rad_half)]
        # histogram verdeelt de radiele distributie in bins, deze bins
        # worden bepaald door ons interval
        fractie = list(histogram(baan_rad_half, interval, weights=weights)[0])
        rad_distr_E.append(fractie)
        # plt.plot(tuple(histogram(baan_rad_half, interval, weights=weights)[1][:-1]), fractie)
        # bovenstaande plot kan getest worden voor 1 E-L koppel
        # run deze plot aub niet in iteratie met L'en en E's, dat gaat
        # poepeloeri lang duren
    # nu is rad_distr_E een lijst met als elk element de radiele distributie
    # (op hun beurt tuples) bij 1 L-waarde. Elk element van rad_distr_E
    # heeft dus 20 elementen (We verdelen in 20 L-waarden) en elk element
    # is een tuple van de radiele verdeling bij een E-L koppel
    return rad_distr_E


def rad_distr_tot(mass_frac, i):
    r_max = r_mass(mass_frac)
    # i is het aantal delen dat we de r_max opdelen
    # een interval opgesteld van 0 tot r_max in 100 stukjes
    sterren_fractie = list()
    # nu gaan we rad_distr_e sommeren ove e om zo een verdeling van sterren
    # te krijgen
    for e in linspace(0.4, 0.95, 40):
        distr = list(rad_distr_e(r_max, e, i))
        sterren_fractie.append(distr)
    return sterren_fractie
    # lijst met elementen [E], met E = [L] en L = een tuple van de radiele
    # distibutie bij 1 E-L koppel. E is dus een lijst van verschillende L'en
    # en sterren_fractie van verschillende E's

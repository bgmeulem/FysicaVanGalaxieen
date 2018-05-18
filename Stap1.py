import numpy
from scipy.integrate import odeint
from scipy.optimize import fmin
from scipy.optimize import brentq
import matplotlib.pyplot as plt

from scipy import interpolate


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
    return -1/((r+1)**2)


def BindPotDer2(r):
    return 2/((r+1)**2)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# tweede stap: Routines om heen en weer te gaan tussen
# Energie E en draaimoment L, pericentrumafstand en
# apocentrumafstand


def aperi(E, L):
    aperi_list = []

    def f(r):
        return(2*E*(r**3) + (2*E - 2)*(r**2) + r*L**2 + L**2)
    if L < 10**(-4):
        # L is praktisch 0: ster oscilleert of zit stil
        aperi_list.append(float(0.0))
        if E != 0:  # ster oscilleert
            aperi_list.append((1-E)/E)
        else:  # ster zit stil
            aperi_list.append(float(0.0))
    else:  # L != 0: cirkelbaan of ellips
        minimum = fmin(f, 10**(-5), xtol=0.000001, disp=False)[0]
        if L == findL(E):
            # cirkelbaan
            aperi_list.append(minimum)
            return aperi_list
        # geen cirkelbaan : 2 nulpunten
        else:
            aperi_list.append(brentq(f, 10**(-6), minimum))
            aperi_list.append(brentq(f, minimum, r_mass(0.99)))
    return aperi_list

# Volgende functies geven minimum en maximum oplossing van 1.112


def aphelium(E, L):
    aperium = list(aperi(E, L))
    return float(max(aperium))


def perihelium(E, L):
    aperium = list(aperi(E, L))
    return float(min(aperium))


def energie(ap, peri):
    return((ap + peri + ap*peri)/((ap + peri)*(ap + 1)*(peri+1)))


def draaimoment(ap, peri):
    if ap == 0:
        return 0
    return numpy.sqrt((2*(ap**2)*(peri**2))/((ap + peri)*(ap + 1)*(peri + 1)))

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
        return 2 / (2*(-E + BindPot(r)) - (L**2 / r**2))**0.5

    if peri != apo:
        periode = abs(integrate.quad(functie, peri, apo)[0])
        if peri == 0:
            return 2*periode
        else:
            return periode
    elif peri == apo:
        return (numpy.pi*2*(apo**2))/L

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Stap 3: Rosettebanen integreren


def BaanInt(apo, peri, stapjes=80):
    L = draaimoment(apo, peri)
    # stapjes moet een even getal zijn want bij radiele oscillatie wordt
    # dit in 2 gesplitst
    stapjes = (stapjes//2)*2
    # Er zijn 4 gevallen waarin onderscheid moet gemaakt worden nl. cirkelbanen,
    # radiele oscillaties, sterren stil staand in het centrum en dan nog de banen
    # die zowel rond het centrum gaan als radiele bewegen

    #  ster staat stil in het centrum: peri=apo=o, L=0
    if peri == apo == 0:
        # De ster staat stil en beweegt niet en heeft een periode die oneindig
        # groot is
        f0 = [0, 0, 0]
        t = []
        oplossingen = []
        for k in range(0, stapjes):
            t.append(k)
            oplossingen.append(f0)

    # De ster beweegt op cirkelbanen: peri=apo>0 en L>0, er moet niet meer gekeken
    # worden of deze niet 0 zijn, dit is hiervoor reeds gebeurd
    elif peri == apo:
        # Men zal nu niet meer de radiele periode moeten volgen maar circulaire
        # periode
        # Hier zal enkel de hoek veranderen
        periode = T_rad(apo, peri)
        f0 = [apo + 10**(-3), 0, 0]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [0, L/(r**2), 0]
            return dfdt
        t = numpy.linspace(0, periode, stapjes)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))

    # De standaard banen: 0<peri, 0<apo en L>0
    elif L > 0:
        # De ster begint in zijn apohelium met hoek=0 en dat is een keerpunt
        # van de snelheid dus is v_r = 0
        f0 = [apo + 10**(-3), 0, 0]
        # integratie werkt niet op apo zelf, v*t = 0

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, L/(r**2), (((L**2)/(r**3)) + BindPotDer1(r))]
            return dfdt
        # Men zal een periode bekijken in gelijke stapjes
        t = numpy.linspace(0, T_rad(peri, apo), stapjes)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))

    # De laatste beweging wordt opgesplitst in 2 delen nl. van apo naar centrum
    # en van centrum terug naar apo maar aan de andere kant van het centrum
    else:
        # Deel 1: van apo naar centrum
        periode = T_rad(apo, peri)/2
        f0 = [apo + 10**(-3), 0, 0]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, 0, BindPotDer1(r)]
            return dfdt
        t1 = numpy.linspace(0,  periode, stapjes//2)
        oplossingen1 = odeint(baanvergelijkingen, f0, t1, args=(L,))

        # Deel 2: van centrum naar apo
        f0 = [0, numpy.pi, -oplossingen1[-1][2]]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, 0, BindPotDer1(r)]
            return dfdt
        t2 = numpy.linspace(periode, 2*periode, stapjes//2)
        oplossingen2 = odeint(baanvergelijkingen, f0, t2, args=(L,))

        # Combinatie van de 2 delen
        oplossingen = oplossingen1.tolist() + oplossingen2.tolist()
        t = t1.tolist() + t2.tolist()

    # De gevonden oplossingen en de time stamp worden nu elk in afzonderlijke
    # sublijsten gestoken van 1 mainlijst
    # Deze ziet er als volgt uit [[tijd][radius][hoek][radiele snelheid]]

    tijd, radius, hoek, radiele_snelheid = [], [], [], []
    for element in range(0, stapjes):
        tijd.append(t[element])
        radius.append(oplossingen[element][0])
        hoek.append(oplossingen[element][1])
        radiele_snelheid.append(oplossingen[element][2])
    return [tijd, radius, hoek, radiele_snelheid]

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Stap 4: voor verschillende r'en de E en L berekenen
# Dan een fit maken zodat men voor willekeurige E de L weet

# een functie die de r geeft bij een gegeven massa, uitgedrukt in decimalen
# Totale massa is 1
# Dus bij 99% van de totale massa is n = 0.99


def r_mass(n):
    r_half = 1 + 2**0.5
    r_max = (((2*n)**0.5)*r_half) / (1 + (1 - (2*n)**0.5)*r_half)
    # bij M = 0.99 is dit ongeveer 198.5
    return r_max


def ListELcouples(r_max):
    r_max = r_mass(0.99)
    E = []
    L = []
    for r in numpy.linspace(r_max/1000, r_max, 1000):
        orb_mom = draaimoment(r, r)
        e = energie(r, r)
        E.insert(0, e)
        L.insert(0, orb_mom)
        # E moet in stijgende volgorde zijn voor de komende plotfunctie
        # E daalt bij hogere r, vandaar insert
    return [E, L]
    # returnt een lijst met als eerste element een lijst van alle E waarden
    # en als tweede de L-waarden


# print(ListELcouples(200))
# print(draaimoment(100))
# print(energie(100))


def findL(E):
    r_max = r_mass(0.99)  # 99% van de totale massa wordt beschouwd
    couples = ListELcouples(r_max)
    E_list = couples[0]
    L_list = couples[1]
    spl = interpolate.UnivariateSpline(E_list, L_list, s=0)
    # plt.plot(numpy.arange(0, 1, 0.0001), spl(numpy.arange(0, 1, 0.0001)))
    # plt.plot(E_list, L_list)
    # plot neemt kwadratisch af naar 0, zoals het hoort
    return spl.__call__(E)  # returnt de waarde van de fit op positie E

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Stap 5: interval 0:r_max opdelen in gelijke intervallen
# Massatoename berekenen
# Integraalruimte opdelen in grid (E_k, L_k)


def interval_r(r_max, n):
    r_part = []
    i = 0
    while i <= r_max:
        r_part.append(i)
        i += r_max/n
    return r_part


def mass_increase(r1, r2):
    if r1:
        return mass(r2) - mass(r1)
    else:
        return mass(r1)


def rad_distr_e(r_max, e, i=100):
    # een radiele distributie voor 1 zekere e
    # het straal-interval wordt standaard verdeeld in 100 stukjes
    interval = numpy.linspace(0, r_mass(0.99), i)
    rad_distr_E = []
    for l in numpy.linspace(0, findL(e) - 10**(-2), 20):
        # Draaimoment bij cirkelbaan is steeds de maximale voor een
        # bepaalde energie
        apo = aphelium(e, l)
        peri = perihelium(e, l)
        baan_rad = BaanInt(apo, peri)[1]
        baan_rad_half = baan_rad[:len(baan_rad)/2]
        # histogram verdeelt de radiële distributie in bins, deze bins
        # worden bepaald door ons interval
        fractie = numpy.histogram(baan_rad_half, bins=interval)[0]
        # deze telt gewoon hoeveel stralen er in een bepaald r_interval
        # zitten, dit dient nog genormeerd te worden zodat de som van
        # alle bins de periode geeft:
        fractie_norm = [x/len(baan_rad_half) for x in fractie]
        rad_distr_E.append(fractie_norm)
    # elke ster krijgt een lijst met de fractie van tijd dat ze
    # doorbrengt in een bepaald r-interval (in i stukken opgedeeld)


def rad_distr_tot(r_max, i=100):
    # i is het aantal delen dat we de r_max opdelen
    # een interval opgesteld van 0 tot r_max in 100 stukjes
    rad_distr_tot = []
    for e in numpy.linspace(10**(-1), 0.9, 20):
        rad_distr_tot.append(rad_distr_e(r_max, e, i))
    return rad_distr_tot

# print(aphelium(0.1, 0.00001))
# print(BaanInt(aphelium(0.1, findL(0.1)), (0.1, findL(0.1)))[1])


plt.plot(ListELcouples(r_mass(0.99))[0], ListELcouples(r_mass(0.99))[1])
print(findL(0.5))
apo = aphelium(0.5, 0.1)
peri = perihelium(0.5, 0.1)
baan_rad = BaanInt(apo, peri)[1]
print(baan_rad)
print(findL(0.1))
distr = rad_distr_tot(r_mass(0.99))
for element in distr:
    print(element)

# t, radius, hoek, snelheid = BaanInt(0, 0)
# plt.plot(t, radius, 'b', label='radius(t)')
# plt.plot(t, hoek, 'g', label='angle(t)')
# plt.plot(t, snelheid, 'r', label='radial velocity(t)')
# plt.legend(loc='best')
# plt.xlabel('t')
# plt.grid()
# plt.show()

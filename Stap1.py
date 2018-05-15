import numpy
from scipy.integrate import odeint
from scipy.optimize import fmin, brentq
import matplotlib.pyplot as plt
from scipy import interpolate
# eerste stap: grootheden


# kleine range functie voor floats
def frange(x, y, jump):
    while x < y:
        yield x
        x += jump


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


# tweede stap: Routines om heen en weer te gaan tussen
# Energie E en draaimoment L, pericentrumafstand en
# apocentrumafstand


def aperi(E, L):
    def f(r):
        return ((L**2)/(2*r**2) - BindPot(r) + E)
    minimum = fmin(f, 10)
    print(minimum[0])
    aperi = []
    aperi.append(scipy.optimize.brentq(f, 0.01, minimum[0]))
    if 10**(-3) < minimum < 10**3:
        # cirkelbaan
        return aperi
    aperi.append(scipy.optimize.brentq(f, minimum[0], r_mass(0.99)))
    return aperi
# Volgende functies geven minimum en maximum oplossing van 1.112


def aphelium(E, L):
    aperium = list(aperi(E, L))
    return max(aperium)


def perihelium(E, L):
    aperium = list(aperi(E, L))
    return min(aperium)


def energie(ap, peri):
    a = numpy.array([[2*(ap**3)+2*(ap**2), ap+1], [2*(peri**3)+2*(peri**2), peri + 1]])
    b = numpy.array([2*(ap**2), 2*(peri**2)])
    return(numpy.linalg.solve(a, b)[0])


def draaimoment(ap, peri):
    a = numpy.array([[2*(ap**3)+2*(ap**2), ap+1], [2*(peri**3)+2*(peri**2), peri + 1]])
    b = numpy.array([2*(ap**2), 2*(peri**2)])
    return(numpy.sqrt(numpy.linalg.solve(a, b)[1]))

# Radiele periode volledig bepaald door peri-en apoheleum
# want men kan via peri en apo de energie en draaimoment bepalen


def T_rad(apo, peri):
    import scipy.integrate as integrate
    E = energie(apo, peri)
    L = draaimoment(apo, peri)
# integrate.quad neemt enkel functie objecten of methoden aan

    def functie(r):
        return 2 / (2*(-E + BindPot(r)) - (L**2 / r**2))**0.5
    return abs(integrate.quad(functie, peri, apo)[0])

# Stap 3: Rosettebanen integreren
# waiting for Triss to finish


def BaanInt(apo, peri):
    # maak onderscheid tussen een situatie bestaande uit
    # enkel radiele oscillaties en werkelijke banen rondom het center
    L = draaimoment(apo, peri)
    if L > 10**(-5):
        # L is niet 0
        # introduceer de beginvoorwaarden
        f0 = [apo, 0, 0]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, L/(r**2), (((L**2)/(r**3)) + BindPotDer1(r))]
            return dfdt

        # men zal kijken waar de ster zich op zijn baan bevivindt gedurende 1 periode T
        # met als tijdstapjes T/80
        t = numpy.linspace(0, T_rad(peri, apo), 81)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))

    # nu hetzelfde maar in het geval van L = 0
    else:
        # in dit geval geldt onze formule voor T_rad niet, eerst bepaalt men
        # deze dus
        # nu hetzelfde verhaal als hierboven
        def baanvergelijkingen(f, t, E, L):
            r, phi, v_r = f
            dfdt = [v_r, 0, BindPotDer1(r)]
            return dfdt
        t = numpy.linspace(0,  T_rad(apo, peri), 81)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))
    return oplossingen

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


def ListELcouples(r_max, r_begin=0):
    E = []
    L = []
    for r in frange(r_begin + (r_max - r_begin)/1000, r_max + (r_max - r_begin)/1000, (r_max - r_begin)/1000):
        if r_begin:
            M = mass(r_max) - mass(r_begin)
        else:
            M = mass(r_max)
        orb_mom = (M*r)**(0.5)
        e = BindPot(r) + M/(2*r)
        E.insert(0, e)
        L.insert(0, orb_mom)
        # E moet in stijgende volgorde zijn voor de komende plotfunctie
        # E daalt bij hogere r, vandaar insert
    return [E, L]
    # returnt een lijst met als eerste element een lijst van alle E waarden
    # en als tweede de L-waarden


def findL(E):
    r_max = r_mass(0.99)  # 99% van de totale massa wordt beschouwd
    couples = ListELcouples(r_max)
    E_list = couples[0]
    L_list = couples[1]
    spl = interpolate.UnivariateSpline(E_list, L_list)
    # plt.plot(E_list, L_list)
    # plot neemt kwadratisch af naar 0, zoals het hoort
    return spl.__call__(E)  # returnt de waarde van de fit op positie E

# Stap 5: interval 0:r_max opdelen in gelijke intervallen
# Massatoename berekenen
# Integraalruimte opdelen in grid (E_k, L_k)


def interval_r(r_max, n):
    r_part = []
    i = 0
    while i <= r_max:
        r_part.append(n)
        i += r_max/n
    return r_part


def mass_increase(r1, r2):
    if r1:
        return mass(r2) - mass(r1)
    else:
        return mass(r1)


def rad_distr(r_max, i=100):
    i = 100
    # i is het aantal delen dat we de r_max opdelen
    # Energie gaat van 1 (r = 0) naar 0 (r = oneindig)
    # delen we deze op in 50 stukjes:
    interval = interval_r(r_mass(0.99), i)
    # een interval opgesteld van 0 tot r_max in 100 stukjes
    sterren_fractie = []
    for e in frange(0.05, 1.05, 0.05):
        # E gaat van 0 naar 1, delen we op in stapjes van 20
        # nu komt voor elk 1/20e van de Energie een interval overeen van
        # mogelijke L'en. Deze delen we ook op in 20 stukjes, genoemd m:
        L_begin = findL(e)
        L_eind = findL(e + 1/20 - 1/400)
        # we hebben de L nodig vlak voor de volgende e (e + 1/20)
        for l in frange(L_begin, L_eind + L_eind/20, (L_begin - L_eind)/20):
            # eerste e is altijd een cirkelbaan, omdat de functie findL
            # werkt met een fit voor cirkelbanen. L verhoogt hierna
            # dus geen cirkelbanen meer
            apo = aphelium(e, l)
            peri = perihelium(e, l)
            baan_rad = BaanInt(apo, peri)[0]
            baan_rad_half = baan_rad[:len(baan_rad)/2]
            # volgende code is erg afhankelijk van de resolutie van
            # baan_int en het interval
            # eerst zoeken wat de laagste waarde is van r in baan_rad_half
            # alvorens we vergelijken met interval
            k = 0
            while baan_rad_half[0] < interval[k]:
                k += 1
            # we gaan ervan uit dat baan_rad_half geen element heeft groter dan
            # r_max
            # nu is baan_rad_half[0] sowieso minstens gelijk aan interval[i]
            for element in baan_rad_half:
                fractie = i*[0]
                # elke ster krijgt een lijst met de fractie van tijd dat ze
                # doorbrengt in een bepaald r-interval (in i stukken opgedeeld)
                if element < interval[k+1]:
                    print(T_rad(apo, peri))
                    fractie[k] += T_rad(apo, peri)/81
                    print(fractie[k])
                else:
                    # we gaan ervan uit dat de resolutie van baan_rad_half
                    # kleiner is dan die van interval, dus als baan_rad_half[3]
                    # groter is dan interval[i], is die sowieso kleiner dan
                    # interval[i+1]
                    i += 1
                    fractie[i] += T_rad(apo, peri)/81
                sterren_fractie.append(fractie)
                # deze lijst wordt toegevoegd aan de totale lijst genaamd
                # sterren_fractie. Elke 20 waarden komen overeen met 1 E-waard
                # en 20 verschillende L-waarden. Deze lijst zou 400 elementen
                # moeten bevatten
    print(sterren_fractie)
    return sterren_fractie



print(energie(1.5, 0.5))
print(draaimoment(1.5, 0.5))
print(aphelium(0.366666666667, 0.387298334621))
print(perihelium(0.366666666667, 0.387298334621))
# print(T_rad(1.5, 0.5))
test = BaanInt(1.5, 0.5)
# print(test)
t = numpy.linspace(0, 9.145870544399036, 81)

plt.plot(t, test[:, 0], 'b', label='radius(t)')
plt.plot(t, test[:, 1], 'g', label='angle(t)')
plt.plot(t, test[:, 2], 'r', label='radial velocity(t)')

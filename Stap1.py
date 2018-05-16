import numpy
from scipy.integrate import odeint
from scipy.optimize import fmin
from scipy.optimize import brentq
#import matplotlib.pyplot as plt

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
    minimum = fmin(f, 0.1)
    aperi = []
    aperi.append(brentq(f, 0.01, minimum[0]))
    if -10**(-3) < f(minimum[0]) < 10**(-3):
        # cirkelbaan
        aperi.append(minimum[0])
        return aperi
    aperi.append(brentq(f, minimum[0], r_mass(0.99)))
    return aperi

# Volgende functies geven minimum en maximum oplossing van 1.112

def aphelium(E, L):
    aperium = list(aperi(E, L))
    return max(aperium)


def perihelium(E, L):
    aperium = list(aperi(E, L))
    return min(aperium)


def energie(ap, peri=0):
    if not peri:
        E = (BindPot(ap) - mass(ap)/(2*ap))
        return E
    else:
        a = numpy.array([[2*(ap**3)+2*(ap**2), ap+1], [2*(peri**3)+2*(peri**2), peri + 1]])
        b = numpy.array([2*(ap**2), 2*(peri**2)])
        return(numpy.linalg.solve(a, b)[0])
    


def draaimoment(ap, peri=0):
    if not peri:
        return ((mass(ap)*ap)**(0.5))
    else:
        a = numpy.array([[2*(ap**3)+2*(ap**2), ap+1], [2*(peri**3)+2*(peri**2), peri + 1]])
        b = numpy.array([2*(ap**2), 2*(peri**2)])
        return(numpy.sqrt(numpy.linalg.solve(a, b)[1]))

# Radiele periode volledig bepaald door peri-en apoheleum
# want men kan via peri en apo de energie en draaimoment bepalen

#in deze functie zit ook nog de periode voor een cirkelbaan, dit om de dingen wat te bundelen
def T_rad(apo, peri):
    import scipy.integrate as integrate
    E = energie(apo, peri)
    L = draaimoment(apo, peri)
# integrate.quad neemt enkel functie objecten of methoden aan

    def functie(r):
        return 2 / (2*(-E + BindPot(r)) - (L**2 / r**2))**0.5
    return abs(integrate.quad(functie, peri, apo)[0]) if peri != apo != 0 else numpy.pi*2*apo**2/L 

# Stap 3: Rosettebanen integreren

def BaanInt(apo, peri, stapjes=80):
    L = draaimoment(apo, peri)
#stapjes moet een even getal zijn want bij radiele oscillatie wordt dit in 2 gesplitst
    stapjes = (stapjes//2)*2
#Er zijn 4 gevallen waarin onderscheid moet gemaakt worden nl. cirkelbanen, radiele oscillaties,
#sterren stil staand in het centrum en dan nog de banen die zowel rond het centrum gaan als radiele bewegen

#De standaard banen: 0<peri, 0<apo en L>0
    if L > 10**(-5):    
#De ster begint in zijn apohelium met hoek=0 en dat is een keerpunt van de snelheid dus is v_r = 0
        f0 = [apo, 0, 0]
        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, L/(r**2), (((L**2)/(r**3)) + BindPotDer1(r))]
            return dfdt
#Men zal een periode bekijken in gelijke stapjes
        t = numpy.linspace(0, T_rad(peri, apo), stapjes)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))

#De ster staat stil in het centrum: peri=apo=o, L=0  
    elif peri == apo == 0:
#De ster staat stil en beweegt niet en heeft een periode die oneindig groot is
        f0 = [0, 0, 0]
        t = []
        oplossingen = []
        for k in range(0,stapjes):
            t.append(k)
            oplossingen.append(f0)

#De ster beweegt op cirkelbanen: peri=apo>0 en L > 0, er moet niet meer gekeken worden of deze niet 0 zijn, dit is hiervoor reeds gebeurt
    elif peri == apo:
#Men zal nu niet meer de radiele periode moeten volgen maar circulaire periode
#Hier zal enkel de hoek veranderen
        periode = T_rad(apo, peri)
        f0 = [apo, 0, 0]
        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [0, L/(r**2), 0]
            return dfdt
        t = numpy.linspace(0, periode, stapjes)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))
        
#De laatste beweging wordt opgesplitst in 2 delen nl. van apo naar centrum
#en van centrum terug naar apo maar aan de andere kant van het centrum
    else:
#Deel 1: van apo naar centrum
        periode = T_rad(apo, peri)/2
        f0 = [apo, 0, 0]
        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, 0, BindPotDer1(r)]
            return dfdt
        t1 = numpy.linspace(0,  periode, stapjes//2)
        oplossingen1 = odeint(baanvergelijkingen, f0, t1, args=(L,))
        
#Deel 2: van centrum naar apo        
        f0 = [0, numpy.pi, -oplossingen1[-1][2]]
        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, 0, BindPotDer1(r)]
            return dfdt
        t2 = numpy.linspace(periode, 2*periode, stapjes//2)
        oplossingen2 = odeint(baanvergelijkingen, f0, t2, args=(L,))

#Combinatie van de 2 delen
        oplossingen = oplossingen1.tolist() + oplossingen2.tolist()
        t = t1.tolist() + t2.tolist()

#De gevonden oplossingen en de time stamp worden nu elk in afzonderlijke sublijsten gestoken van 1 mainlijst
#Deze ziet er als volgt uit [[tijd][radius][hoek][radiele snelheid]]

    tijd, radius, hoek, radiele_snelheid = [], [], [], []
    for element in range(0, stapjes):
        tijd.append(t[element])
        radius.append(oplossingen[element][0])
        hoek.append(oplossingen[element][1])
        radiele_snelheid.append(oplossingen[element][2])
    return [tijd, radius, hoek, radiele_snelheid]

    
 
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
    E = []
    L = []
    for r in numpy.linspace(r_max/1000, r_max, 1000):
        orb_mom = draaimoment(r)
        e = energie(r)
        E.insert(0, e)
        L.insert(0, orb_mom)
        # E moet in stijgende volgorde zijn voor de komende plotfunctie
        # E daalt bij hogere r, vandaar insert
    return [E, L]
    # returnt een lijst met als eerste element een lijst van alle E waarden
    # en als tweede de L-waarden
print(ListELcouples(200))
print(draaimoment(100))
print(energie(100))

def findL(E):
    r_max = r_mass(0.99)  # 99% van de totale massa wordt beschouwd
    couples = ListELcouples(r_max)
    E_list = couples[0]
    L_list = couples[1]
    spl = interpolate.UnivariateSpline(E_list, L_list)
    # plt.plot(E_list, L_list)
    # plot neemt kwadratisch af naar 0, zoals het hoort
    return spl.__call__(E)  # returnt de waarde van de fit op positie E
print(findL(energie(100)))
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


def rad_distr(r_max, i=100):
    i = 100
    # i is het aantal delen dat we de r_max opdelen
    interval = interval_r(r_mass(0.99), i)
    # een interval opgesteld van 0 tot r_max in 100 stukjes
    sterren_fractie = []
    for e in numpy.linspace(0,1,20):
        # E gaat van 0 naar 1, delen we op in stapjes van 20
        # we hebben de L nodig vlak voor de volgende e (e + 1/20)
        for l in frange(0, findL(e), 20):
            # Draaimoment bij cirkelbaan is steeds de maximale voor een bepaalde energie
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
                # sterren_fractie. Elke 20 waarden komen overeen met 1 E-waarde
                # en 20 verschillende L-waarden. Deze lijst zou 400 elementen
                # moeten bevatten
    return sterren_fractie

#t, radius, hoek, snelheid = BaanInt(0, 0)
#plt.plot(t, radius, 'b', label='radius(t)')
#plt.plot(t, hoek, 'g', label='angle(t)')
#plt.plot(t, snelheid, 'r', label='radial velocity(t)')
#plt.legend(loc='best')
#plt.xlabel('t')
#plt.grid()
#plt.show()
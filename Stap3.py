# Stap 3: Rosettebanen integreren
from Stap2 import draaimoment
from Stap1 import BindPotDer1
from numpy import linspace
from scipy.integrate import odeint
from Stap2 import T_rad
from numpy import pi
import matplotlib.pyplot as plt


def BaanInt(apo, peri, stapjes=1000):
    L = draaimoment(apo, peri)
    # stapjes moet een even getal zijn want bij radiele oscillatie wordt
    # dit in 2 gesplitst
    stapjes = (stapjes//2)*2
    # Er zijn 4 gevallen waarin onderscheid moet gemaakt worden nl. cirkelbanen
    # radiele oscillaties, sterren stilstaand in het centrum en dan nog de
    # banen die zowel rond het centrum gaan als radiele bewegen
    # De standaard banen: 0<peri, 0<apo en L>0
    if 0 < peri and peri != apo:
        # De ster begint in zijn apohelium met hoek=0 en dat is een keerpunt
        # van de snelheid dus is v_r = 0
        f0 = [apo, 0, 0]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f

            dfdt = [v_r, L/(r**2), (((L**2)/(r**3)) + BindPotDer1(r))]
            return dfdt
        # Men zal een periode bekijken in gelijke stapjes
        t = linspace(0, T_rad(apo, peri), stapjes)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))
    #  ster staat stil in het centrum: peri=apo=o, L=0
    elif peri == apo == 0:
        # De ster staat stil en beweegt niet en heeft een periode die oneindig
        # groot is
        f0 = [0, 0, 0]
        t = []
        oplossingen = []
        for k in range(0, stapjes):
            t.append(k)
            oplossingen.append(f0)

    # De ster beweegt op cirkelbanen: peri=apo>0 en L>0, er moet niet meer
    # gekeken worden of deze niet 0 zijn, dit is hiervoor reeds gebeurd
    elif peri == apo:
        # Men zal nu niet meer de radiele periode moeten volgen maar circulaire
        # periode
        # Hier zal enkel de hoek veranderen
        periode = T_rad(apo, peri)
        f0 = [apo, 0, 0]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [0, L/(r**2), 0]
            return dfdt
        t = linspace(0, periode, stapjes)
        oplossingen = odeint(baanvergelijkingen, f0, t, args=(L,))

    # De laatste beweging wordt opgesplitst in 2 delen nl. van apo naar centrum
    # en van centrum terug naar apo maar aan de andere kant van het centrum
    else:
        # Deel 1: van apo naar centrum
        periode = T_rad(apo, peri)/2
        f0 = [apo, 0, -10**(-3)]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, 0, BindPotDer1(r)]
            return dfdt
        t1 = linspace(0,  periode, stapjes//2)
        oplossingen1 = odeint(baanvergelijkingen, f0, t1, args=(L,))

        # Deel 2: van centrum naar apo
        f0 = [0, pi, -oplossingen1[-1][2]]

        def baanvergelijkingen(f, t, L):
            r, phi, v_r = f
            dfdt = [v_r, 0, BindPotDer1(r)]
            return dfdt
        t2 = linspace(periode, 2*periode, stapjes//2)
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


test = BaanInt(1.5, 0.5, 1000)
plt.plot(test[0], test[1], 'b', label='radius(t)')
plt.plot(test[0], test[2], 'g', label='angle(t)')
plt.plot(test[0], test[3], 'r', label='radial velocity(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

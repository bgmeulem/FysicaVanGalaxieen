# import matplotlib.pyplot as plt


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


# t, radius, hoek, snelheid = BaanInt(0, 0)
# plt.plot(t, radius, 'b', label='radius(t)')
# plt.plot(t, hoek, 'g', label='angle(t)')
# plt.plot(t, snelheid, 'r', label='radial velocity(t)')
# plt.legend(loc='best')
# plt.xlabel('t')
# plt.grid()
# plt.show()

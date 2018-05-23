import matplotlib.pyplot as plt
from Stap6 import m_k
from Stap5 import mass_increase
from Stap5 import rad_distr_e
from Stap4 import r_mass
from numpy import linspace
from Stap6 import BijdrageL
import time
def verwachte_massaverdeling (massfrac, stapjes = 100):
    massa_spreiding = mass_increase(massfrac, stapjes)
    r_max = r_mass(massfrac)
    straal_interval = linspace(0,r_max,stapjes).tolist()
    return straal_interval, massa_spreiding, straal_interval[1]

def verkregen_massaverdeling (massfrac, stapjes=100):
    #Bepalen van de r_i bins
    r_max = r_mass(massfrac)
    straal_interval = linspace(0,r_max,stapjes).tolist()
    
    #Bepalen van de sum over k van (m_k*f_i(E_k))
    #dit is een lijst met als element de m_k horende bij de E_k
    m_k_list = m_k(massfrac, stapjes)
    #De lijst van de E_k's
    E_k_list = linspace(0.0001, 0.999, stapjes).tolist()
    #Opstellen van de lijst met als elementen m_k*f_i(E_k)
    #Met andere woorden de elementen zijn de verkregen massa 
    #aan sterren horende bij energie E_K op straal r_i
    mass_op_straal_list = [0 for x in range(len(straal_interval))]
    #In essentie de som over k
    for i in range(len(straal_interval)):
        m_k_float = float(m_k_list[i])
        #Dit geeft ons de f_i's horende bij een E_k
        r_verdeling_van_E_k = BijdrageL(E_k_list[i], massfrac, stapjes)
        #Nu wordt mass_op_straal_lijst aangepast door voor ieder straalinterval van deze lijst 
        #de waarde aan te vullen met de fractie van de massa horende bij een bepaalde energie
        #bekijkt in essentie elke straal r_i apart op deze manier
        for element in range(len(straal_interval)):
            mass_op_straal_list[element] += (m_k_float * r_verdeling_van_E_k[element])
    return straal_interval, mass_op_straal_list, straal_interval[1]
 
#Om te plotten dient men gewoon de ''' tekens weg te doen en het programma te runnen
#merk op dat dit even kan duren
'''
m_k_list = m_k(0.99, 300)
E_k_list = linspace(0.0001, 0.999, 300).tolist()
plt.plot(E_k_list, m_k_list, 'b', label='m_k(E_k)')
plt.legend(loc='best')
plt.xlabel('E_k')
plt.grid()
plt.show()
'''
'''
from matplotlib.pyplot import *
#stapjes mag nie te groot worden, dit duurt al gigantisch lang
tussenvar = verkregen_massaverdeling(0.9, stapjes = 100)
for x in range(len(tussenvar)):
    tussenvar[0][x] += tussenvar[2]/2
x,y,breedte = tussenvar[0], tussenvar[1], tussenvar[2]
bar(x,y, width = breedte)
show()
'''

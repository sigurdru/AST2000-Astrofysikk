#DETTE ER EGEN KODE
import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as constants

system = SolarSystem(8913) #lager et solsystem.

class my_solar_system:
    def __init__(self):
        self.eksentrisitet_list = system.eccentricities
        self.perihelvinkel_list = system.aphelion_angles-np.pi  #radianer
        self.planetmasse_list = system.masses                   #masse i solmaess
        self.store_halvakse_list = system.semi_major_axes       #AU
        self.initialvel = np.array(system.initial_velocities)   #initialhastigheter AU/år
        self.initialpos = np.array(system.initial_positions)    #initialposisjoner
        self.G_sol = constants.G_sol                            #gravitasjonskonstant skalert
        self.m_sun = system.star_mass                           #massen til solen i solmasser
        self.m_solar = constants.m_sun                          #solmasse i kg
        self.m_mars = 6.39*10**23                               #massen til marsk i kg
        self.G = constants.G                                    #vanlig gravitasjonskonstant
        self.AU = constants.AU                                  #astronomisk enhet
    def get_info(self):
        """
        Printer masse informasjon om solssystemet.
        """
        print('My system has a {:g} solar mass star with a radius of {:g} kilometers.'
        .format(system.star_mass*1.98892*10**30, system.star_radius))
        for planet_idx in range(system.number_of_planets):
            print('Planet {:d} is a {} planet with a semi-major axis of {:g} AU. Masse: {:g}'
                .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx],system.masses[planet_idx]*1.98892*10**30))
        for planet_idx in range(system.number_of_planets):
            print('planet {:d} omega: {:g} e: {:g}'
                .format(planet_idx, system.aphelion_angles[planet_idx]-np.pi, system.eccentricities[planet_idx]))
        print("initial posisjoner")
        print(system.initial_positions)
        print("initial hastigheter")
        print(system.initial_velocities)
    def analytical_plot(self):
        """
        plotter analytisk plot til planetbanene
        """
        angle_array = np.linspace(0,2*np.pi,1000)
        for i in range(len(self.planetmasse_list)):
            value = self.store_halvakse_list[i]*(1-self.eksentrisitet_list[i]**2)/(1+self.eksentrisitet_list[i]*np.cos(angle_array))
            x_array = value*np.cos(angle_array+self.perihelvinkel_list[i])
            y_array = value*np.sin(angle_array+self.perihelvinkel_list[i])
            plt.plot(x_array,y_array, label="Planet %i" %(i+1))
            plt.plot(self.initialpos[0,i],self.initialpos[1,i], marker='o')
        plt.plot([0],[0], marker='*')
        plt.xlabel("x-posisjon i AU")
        plt.ylabel("y-posisjon i AU")
        plt.title("Analytisk planetbane")
        plt.legend()
        plt.axis("equal")
    def numerical_plot_solsys(self,dt,t_slutt):
        """
        finner numerisk planetbane
        """
        n = int(t_slutt/dt)
        pos_array = np.zeros((2,system.number_of_planets,n))
        vel_array = np.zeros((2,system.number_of_planets,n))
        for j in range(system.number_of_planets):
            pos_array[:,j,0] = self.initialpos[:,j]
            vel_array[:,j,0] = self.initialvel[:,j]
            for i in range(n-1):
                a = -self.G_sol*self.m_sun*(pos_array[:,j,i]/np.linalg.norm(pos_array[:,j,i])**3)
                vel_array[:,j,i+1] = a*dt + vel_array[:,j,i]
                pos_array[:,j,i+1] = vel_array[:,j,i+1]*dt + pos_array[:,j,i]
        return vel_array,pos_array
    def numerical_plot_three(self,dt,n):
        """
        numerisk plot av trelegemeproblemet
        """
        pos_big_star = np.zeros((2,n)) #lager posisjons og hastighetsarrays
        pos_small_star = np.zeros((2,n))
        pos_planet = np.zeros((2,n))

        vel_big_star = np.zeros((2,n))
        vel_small_star = np.zeros((2,n))
        vel_planet = np.zeros((2,n))

        pos_big_star[0,0] = 3*self.AU    #skalerer
        pos_planet[0,0] = -1.5*self.AU  

        vel_planet[1,0] = -1000          #initialposisjoner                  
        vel_small_star[1,0] = 30000
        vel_big_star[1,0] = -7500

        for i in range(n-1):
            big_to_planet = pos_planet[:,i]-pos_big_star[:,i]    #finner vektorene mellom legemene
            small_to_planet = pos_planet[:,i]-pos_small_star[:,i]
            small_to_big = pos_big_star[:,i]-pos_small_star[:,i]

            a_planet = self.G*(4*self.m_solar*(-big_to_planet/np.linalg.norm(big_to_planet)**3)\  
                       + self.m_solar*(-small_to_planet/np.linalg.norm(small_to_planet)**3))        #regner ut akselerasjonen til alle planetene

            a_big_star = self.G*(self.m_mars*(big_to_planet/np.linalg.norm(big_to_planet)**3)\
                       + self.m_solar*(-small_to_big/np.linalg.norm(small_to_big)**3))

            a_small_star = self.G*(self.m_mars*(small_to_planet/np.linalg.norm(small_to_planet)**3)\
                         + 4*self.m_solar*(small_to_big/np.linalg.norm(small_to_big)**3))

            vel_planet[:,i+1] = a_planet*dt + vel_planet[:,i]               #bruker Euler-Cromer på alle planetene
            vel_big_star[:,i+1] = a_big_star*dt + vel_big_star[:,i]
            vel_small_star[:,i+1] = a_small_star*dt + vel_small_star[:,i]

            pos_planet[:,i+1] = vel_planet[:,i+1]*dt + pos_planet[:,i]
            pos_big_star[:,i+1] = vel_big_star[:,i+1]*dt + pos_big_star[:,i]
            pos_small_star[:,i+1] = vel_small_star[:,i+1]*dt + pos_small_star[:,i]
        return pos_planet/self.AU,pos_big_star/self.AU,pos_small_star/self.AU

dt1 = 400
n1 = 10**6
pos_planet, pos_big_star, pos_small_star = my_solar_system().numerical_plot_three(dt1,n1)

times1 = np.linspace(0,int(dt1*n1),n1)

plt.plot(pos_planet[0],pos_planet[1], label="Planet")
plt.plot(pos_big_star[0],pos_big_star[1], label="Big star")
plt.plot(pos_small_star[0],pos_small_star[1], label="Small star")
plt.title("Trelegeme")
plt.xlabel("x-posisjon i AU")
plt.ylabel("y-posisjon i AU")
plt.axis("equal")
plt.legend()
plt.show()

#lager videofil
system.generate_binary_star_orbit_video(times1,pos_planet,pos_big_star,pos_small_star, filename='binary_orbit_video.xml')


#plotter analytisk plot
my_solar_system().analytical_plot()
plt.show()


t_slutt = 15
dt2 = 0.01
pos_array = my_solar_system().numerical_plot_solsys(dt2,t_slutt)[1]
times2 = np.linspace(0,t_slutt,int(t_slutt/dt2))
system.generate_orbit_video(times2, pos_array, filename='orbit_video.xml')


for i in range(system.number_of_planets):
    plt.plot(pos_array[0,i],pos_array[1,i], label="Planet %i" %int(i+1))
    plt.plot(system.initial_positions[0,i], system.initial_positions[1,i], marker='o')    

plt.plot([0],[0], marker='*')
plt.title("Simulert planetbane")
plt.xlabel("x-posisjon i AU")
plt.ylabel("y-posisjon i AU")
plt.axis("equal")
plt.legend()
plt.show()


#finner informasjonen til planetene
my_solar_system().get_info()
"""
Resultat:
My system has a 1.08239 solar mass star with a radius of 751861 kilometers.
Planet 0 is a rock planet with a semi-major axis of 0.680887 AU. Masse: 5.41211e+24
Planet 1 is a rock planet with a semi-major axis of 0.960048 AU. Masse: 4.49684e+24
Planet 2 is a gas planet with a semi-major axis of 4.57042 AU. Masse: 1.36195e+26
Planet 3 is a gas planet with a semi-major axis of 6.10486 AU. Masse: 6.77856e+26
Planet 4 is a rock planet with a semi-major axis of 1.51256 AU. Masse: 1.30022e+23
Planet 5 is a gas planet with a semi-major axis of 3.34493 AU. Masse: 2.99985e+27
Planet 6 is a rock planet with a semi-major axis of 0.499398 AU. Masse: 2.49086e+24
Planet 7 is a gas planet with a semi-major axis of 2.32955 AU. Masse: 1.26516e+26
planet 0 omega: -3.14159 e: 0.0273064
planet 1 omega: 1.31124 e: 0.0186564
planet 2 omega: 2.10767 e: 0.0652503
planet 3 omega: 0.1823 e: 0.0384509
planet 4 omega: -0.476925 e: 0.0360803
planet 5 omega: -1.09926 e: 0.0125299
planet 6 omega: 2.83067 e: 0.0698374
planet 7 omega: -1.2164 e: 0.0072167
initial posisjoner
[[ 0.69947949 -0.5959163  -4.32591143  1.49065548 -1.21320213 -1.56327789 -0.22589133  0.46007707]
 [ 0.          0.7388836   0.62861253  5.81080199 -0.9386583  -2.9295989  -0.43624904 -2.26671946]]
initial hastigheter
[[ 0.         -5.31428591 -0.61245741 -2.58303629  3.3427196   3.19348444   8.03624823  4.22639155]
 [ 7.70855232 -4.15703091 -3.13463824  0.7580127  -4.03607685 -1.66245446  -4.88032314  0.86267508]]
"""



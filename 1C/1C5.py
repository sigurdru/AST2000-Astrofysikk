#DETTE ER EGEN KODE
import numpy as np
import matplotlib.pyplot as plt
import random as random

import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as constants

system = SolarSystem(8913)                                                      #lager et solsystem.

class orbit:
    def __init__(self):
        self.G_sol = constants.G_sol                                            #gravitasjonskonstanten
        self.mass_array = system.masses                                         #array med massene til planetene i solmasser [sol]
        self.m_star = system.star_mass                                                          #massen til stjernen i solmasser [sol]
        self.initial_positions = np.array([np.array(i) for i in system.initial_positions])      #initial posisjonen til planetene [AU]
        self.initial_velocities = np.array([np.array(i) for i in system.initial_velocities])     #initial hastigheter til planetene [AU/år]

    def find_massive_planets(self):
        mass_array_copy = self.mass_array.copy()                                #lager en kopi av mass_array så jeg kan redigere den
        planet_indexes = []
        for _ in range(3):
            index = np.where(mass_array_copy == np.amax(mass_array_copy))       #finner index til største massen
            planet_indexes.append(index)                                        #legger til index
            mass_array_copy[index] = -1                                         #ingen av massene har negativ masse
        return planet_indexes
        
    def numerical_plot(self, n, dt):
        m_star = self.m_star
        planet_indexes = self.find_massive_planets()
        m1 = self.mass_array[planet_indexes[0]][0]                              #når jeg indekserer med planet_indexes får jeg uønsket lister
        m2 = self.mass_array[planet_indexes[1]][0]                              #derfor må jeg fjerne det med [0]
        m3 = self.mass_array[planet_indexes[2]][0]

        v_star = np.zeros((2,n))                                                #definerer hastighets arrays
        v_planet1 = np.zeros((2,n))
        v_planet2 = np.zeros((2,n))
        v_planet3 = np.zeros((2,n))

        v_planet1[0,0] = self.initial_velocities[:,planet_indexes[0]][0,0,0]    #her også får jeg uønsket liter jeg må fjerne
        v_planet1[1,0] = self.initial_velocities[:,planet_indexes[0]][1,0,0]
        v_planet2[0,0] = self.initial_velocities[:,planet_indexes[1]][0,0,0]
        v_planet2[1,0] = self.initial_velocities[:,planet_indexes[1]][1,0,0]
        v_planet3[0,0] = self.initial_velocities[:,planet_indexes[2]][0,0,0]
        v_planet3[1,0] = self.initial_velocities[:,planet_indexes[2]][1,0,0]

        r_star = np.zeros((2,n))                                                #definerer posisjons arrays
        r_planet1 = np.zeros((2,n))
        r_planet2 = np.zeros((2,n))
        r_planet3 = np.zeros((2,n))
        
        r_planet1[0,0] = self.initial_positions[:,planet_indexes[0]][0,0,0]
        r_planet1[1,0] = self.initial_positions[:,planet_indexes[0]][1,0,0]
        r_planet2[0,0] = self.initial_positions[:,planet_indexes[1]][0,0,0]
        r_planet2[1,0] = self.initial_positions[:,planet_indexes[1]][1,0,0]
        r_planet3[0,0] = self.initial_positions[:,planet_indexes[2]][0,0,0]
        r_planet3[1,0] = self.initial_positions[:,planet_indexes[2]][1,0,0]

        for i in range(n-1):

            star_to_planet1 = r_planet1[:,i]-r_star[:,i]                        #finner posisjonsvektorene mellom planetene
            star_to_planet2 = r_planet2[:,i]-r_star[:,i]
            star_to_planet3 = r_planet3[:,i]-r_star[:,i]
            planet1_to_planet2 = r_planet2[:,i]-r_planet1[:,i]
            planet1_to_planet3 = r_planet3[:,i]-r_planet1[:,i]
            planet2_to_planet3 = r_planet3[:,i]-r_planet2[:,i]

            a_star = self.G_sol*(m1*star_to_planet1/np.linalg.norm(star_to_planet1)**3 +\
                                m2*star_to_planet2/np.linalg.norm(star_to_planet2)**3 +\
                                m3*star_to_planet3/np.linalg.norm(star_to_planet3)**3)
            a1 = self.G_sol*(m_star*(-star_to_planet1)/np.linalg.norm(star_to_planet1)**3 +\
                            m2*planet1_to_planet2/np.linalg.norm(planet1_to_planet2)**3 +\
                            m3*planet1_to_planet3/np.linalg.norm(planet1_to_planet3)**3)
            a2 = self.G_sol*(m_star*(-star_to_planet2)/np.linalg.norm(star_to_planet2)**3 +\
                            m1*(-planet1_to_planet2)/np.linalg.norm(planet1_to_planet2)**3 +\
                            m3*(planet2_to_planet3)/np.linalg.norm(planet2_to_planet3)**3)
            a3 = self.G_sol*(m_star*(-star_to_planet3)/np.linalg.norm(star_to_planet3)**3 +\
                            m1*(-planet1_to_planet3)/np.linalg.norm(planet1_to_planet3)**3 +\
                            m2*(-planet2_to_planet3)/np.linalg.norm(planet2_to_planet3)**3)

            v_star[:,i+1] = a_star*dt + v_star[:,i]
            v_planet1[:,i+1] = a1*dt + v_planet1[:,i]
            v_planet2[:,i+1] = a2*dt + v_planet2[:,i] 
            v_planet3[:,i+1] = a3*dt + v_planet3[:,i] 

            r_star[:,i+1] = v_star[:,i+1]*dt + r_star[:,i]
            r_planet1[:,i+1] = v_planet1[:,i+1]*dt + r_planet1[:,i]
            r_planet2[:,i+1] = v_planet2[:,i+1]*dt + r_planet2[:,i]
            r_planet3[:,i+1] = v_planet3[:,i+1]*dt + r_planet3[:,i]
        
        return r_star, r_planet1, r_planet2, r_planet3, v_star

    def observer(self, n, dt, v_drift, inclination):
        star_vel = self.numerical_plot(n,dt)[-1].copy()
        observed = np.sin(inclination)*(star_vel[0]) + v_drift
        t_array = np.linspace(0,int(n*dt),n)
        return t_array, observed

    def add_noise(self, n, dt, v_drift, inclination):
        t_array, observed = self.observer(n,dt,v_drift,inclination)
        v_star = self.numerical_plot(n,dt)[-1]
        v_cm = sum(v_star)/len(v_star)
        v_star = v_star-v_cm
        v_adjust = sum(observed)/len(observed)
        observed = observed-(v_adjust)
        observed_noise = observed.copy()
        mu = 0
        sigma = max(v_star[0])/5
        gauss_array = np.random.normal(mu,sigma,n)
        for i in range(len(observed)):
            noise = gauss_array[i]
            observed_noise[i] += noise
        plt.plot(t_array,observed_noise, label = "inklinasjon = %f" %(inclination))
        plt.plot(t_array,observed, label = "faktisk hastighet")
        

            
t_slutt = 25
dt = 0.01
n = int(t_slutt/dt)
v_drift = 10
inclination_list = [np.pi/2, np.pi/4, np.pi/16, 0]
print(n)

"""
r_star, r_planet1, r_planet2, r_planet3 = orbit().numerical_plot(n,dt)[:-1]

plt.plot(r_star[0],r_star[1], label= "stjerne")
plt.plot(r_planet1[0],r_planet1[1], label = "planet1")
plt.plot(r_planet2[0],r_planet2[1], label = "planet2")
plt.plot(r_planet3[0],r_planet3[1], label = "planet3")
plt.xlabel("posisjon i x-retning[AU]")
plt.ylabel("posisjon i y-retning[AU]")
plt.autoscale("false")
plt.axis("equal")
plt.legend()
plt.show()

for i in inclination_list:
    t_array, observed = orbit().observer(n,dt,v_drift,i)
    plt.plot(t_array,observed, label="inklinasjon: %f" %(i))
plt.xlabel("tid[år]")
plt.ylabel("hastighet[AU/år]")
plt.legend()
plt.show()

for i in inclination_list:
    orbit().add_noise(n,dt,v_drift,i)
    plt.xlabel("tid[år]")
    plt.ylabel("hastighet[AU/år]")
    plt.legend()
    plt.show()
"""
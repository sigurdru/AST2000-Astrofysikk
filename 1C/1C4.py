#DETTE ER EGEN KODE
import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as constants

seed = 8913
system = SolarSystem(seed) #lager et solsystem.


class C4:
    def __init__(self):        
        self.stars = ["star0_4.81.txt", "star1_1.79.txt", "star2_1.02.txt", "star3_1.79.txt", "star4_3.57.txt"] #filnavn
        self.masses = [4.81, 1.79, 1.02, 1.79, 357]  #solmasse        
        self.lambda_0 = 656.28  #[nm]=10^-9m
        self.c = constants.c #[m/s] lysfarten
    def get_values(self):
        """
        Henter ut verdiene fra filene.
        """
        t_list = [] #lager en lste for alle verdiene
        bl_list = []
        flux_list = []
        for i in range(len(self.masses)):
            t_list.append([]) #lager en nestet liste for hver stjerne
            bl_list.append([])
            flux_list.append([])

            infile = open(self.stars[i], "r")
            for line in infile:
                line = line.split()     #leser av filen og fyller alle listene
                t_list[i].append(float(line[0]))
                bl_list[i].append(float(line[1]))
                flux_list[i].append(float(line[2]))
            infile.close()
        t_array = np.array([np.array(i) for i in t_list])   #gjør de til type array
        bl_array = np.array([np.array(i) for i in bl_list])
        flux_array = np.array([np.array(i) for i in flux_list])
        return t_array, bl_array, flux_array
    def find_vel(self):
        t_array, bl_array = C4().get_values()[0:2]  #henter listene jeg trenger

        v_array = self.c*((bl_array-self.lambda_0)/self.lambda_0) #regner ut v_pec
        v_av = np.array([np.mean(v_array[i]) for i in range(len(v_array))]) #regner ut gjennomsnitlig hastighet
        v_egen = v_array[:]-v_av[:]
        return v_egen, t_array
    def plot_light_curves(self,n):
        t_list, bl_array, flux_list = C4().get_values()
        plt.plot(t_list[n],flux_list[n], label ="Observert fluks")
    def plot_closeup(self,n,start, stop):
        t_list, ignore, flux_list = C4().get_values()
        plt.title("stjerne %i" %(n))
        plt.plot(t_list[n][start:stop],flux_list[n][start:stop], label="Mellom indeksering %i og %i" %(start,stop))
        plt.xlabel("tid[dager]")
        plt.ylabel("relativ fluks")
        plt.legend()
    def v_model(self, t, t_0, P, v_r):
        v = v_r*np.cos((2*np.pi/P)*(t-t_0))
        return v
    def minste_kvadraters(self, t0_min, t0_max, P_min, P_max, vr_min, vr_max, stjerne):
        t_array = C4().get_values()[0]
        v_array = C4().find_vel()[0]
        n = 20
        return_list = [[],[],[]]

        for t0min, t0max, Pmin, Pmax, vrmin, vrmax, i in zip(t0_min, t0_max, P_min, P_max, vr_min, vr_max, stjerne):
            t0_array = np.linspace(t0min,t0max,n)
            P_array = np.linspace(Pmin,Pmax,n)
            vr_array = np.linspace(vrmin,vrmax,n)
            dif_list = []
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        v_model_value = self.v_model(t = t_array[i], t_0 = t0_array[j], P = P_array[k], v_r = vr_array[l])
                        dif_list.append(np.sum((v_model_value-v_array[i])**2))
        
            dif_array = np.array(dif_list)
            index = dif_array.argmin()
            t0_index = index//n**2
            P_index = index//n - t0_index*n
            vr_index = index - P_index*n - t0_index*n**2
            return_list[0].append(t0_array[t0_index])
            return_list[1].append(P_array[P_index]) 
            return_list[2].append(vr_array[vr_index])
        return return_list


t0min_list = [3800,1000, 3000]
t0max_list = [4400, 1900, 3750]
Pmin_list = [4000, 3900,4500]
Pmax_list = [5000, 5500, 6000]
stjerne_list = [1, 2, 3]
vrmin_list = [2.5, 7.5, 25]
vrmax_list = [7.5, 12.5, 33]
a = C4().minste_kvadraters(t0_min = t0min_list, t0_max = t0max_list, P_min = Pmin_list, \
                        P_max = Pmax_list, vr_min = vrmin_list, vr_max = vrmax_list, stjerne = stjerne_list)

v_egen, t_array = C4().find_vel()
print(a)

i2 = 0
for i in stjerne_list:
    plt.title("Hastighetsplot til stjerne %i" %(i))
    plt.plot(t_array[i], v_egen[i], label = "Observert hastighet om massesenter")
    plt.plot(t_array[i], C4().v_model(t_array[i],a[0][i2],a[1][i2],a[2][i2]), label="Minste kvadraters metode")
    plt.xlabel("tid[dager]")
    plt.ylabel("hastighet[m/s]")
    plt.legend()
    plt.show()
    i2+=1

lengde2 = 8869          #lengden til flux listen til stjerne 2
lengde3 = 10280         #lengden til flux listen til stjerne 3
C4().plot_closeup(n= 2,start =int(lengde2*0.3),stop= int(lengde2*0.31))
plt.show()
C4().plot_closeup(n=2,start = int(lengde2*0.8),stop=int(lengde2*0.81))
plt.show()
C4().plot_closeup(n=3,start = int(lengde3*0.444),stop=int(lengde3*0.455))
plt.show()
C4().plot_closeup(n=3,start = int(lengde3*0.94),stop=int(lengde3*0.95))
plt.show()



for i in range(5):
    plt.subplot(2,1,1)
    v_egen, t_array = C4().find_vel()
    plt.title("Stjerne %i" %(i))
    plt.plot(t_array[i],v_egen[i], label ="Observert hastighet om massesenter")    
    plt.xlabel("tid[dager]")
    plt.ylabel("hastighet[m/s]")
    plt.legend()
    plt.subplot(2,1,2)
    C4().plot_light_curves(n=i)
    plt.xlabel("tid[dager]")
    plt.ylabel("relativ fluks")
    plt.legend()
    plt.show()


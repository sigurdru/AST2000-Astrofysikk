#DETTE ER EGEN KODE
import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as constants

system = SolarSystem(8913)                                                      #lager et solsystem.

class D6:
    def __init__(self):
        self.days = [0,2,4,6,8,10,13,15,17,19]
        self.c = constants.c
    def get_data(self):
        #Leser av all data og setter det i en stor array
        data_array = np.zeros((len(self.days),1500,2))
        for i in range(len(self.days)):
            data_array[i] = np.loadtxt("spectrum_day%s.txt" %(str(self.days[i])))
        return data_array

    def plot_data(self, x):
        #plotter graf nummer x
        data_array = self.get_data()
        plt.rc('axes', labelsize=12)  
        plt.plot(data_array[x,:,0],data_array[x,:,1], "k", label = "dag %i" %(self.days[x]))
        plt.xlabel("Bølgelengde[nm]")
        plt.ylabel("Flux")
        plt.tight_layout()
        plt.legend()  

    def F_model(self, F_max, F_min, sigma, lambda_center, data_array):
        #regner ut og returnerer F_model
        F_model = F_max + (F_min-F_max)*np.exp(-1*(data_array[:,0]-lambda_center)**2/(2*sigma**2))
        return F_model

    def least_squares(self, F_min_list, sigma_list, lambda_center_list, x):
        #tar in min og max verdier for F, sigma, og lambda, og bruker minste kvadraters på graf x
        data_array = self.get_data()
        F_max = 1
        n = 20
        F_min_array = np.linspace(F_min_list[0], F_min_list[1],n)
        sigma_array = np.linspace(sigma_list[0], sigma_list[1],n)
        lambda_center_array = np.linspace(lambda_center_list[0], lambda_center_list[1],n)

        sum_list = []
        for F_min in F_min_array: 
            for sigma in sigma_array:
                for lambda_center in lambda_center_array: #regner ut summen av minste kvadraters
                    sum_list.append(np.sum((data_array[x,:,1] - self.F_model(F_max, F_min,sigma,lambda_center,data_array[x]))**2))

        index = np.array(sum_list).argmin() #finner indekser
        F_min_index = index//n**2
        sigma_index = index//n - F_min_index*n
        lambda_center_index = index - sigma_index*n - F_min_index*n**2
        return F_min_array[F_min_index], sigma_array[sigma_index], lambda_center_array[lambda_center_index]
    
    def auto_find(self, x):
        data_array = self.get_data()[x]
        #observerer at alle har sirka samme sigma, burde også forvente dette siden det er samme stjerne med samme temperatur
        F_min_index = data_array[:,1].argmin()
        F_min_list = [data_array[F_min_index,1]+0.05,data_array[F_min_index,1]+0.15]
        sigma_list = [0.001,0.02]
        lambda_center_list = [data_array[F_min_index,0]-0.01,data_array[F_min_index,0]+0.01]
        return F_min_list, sigma_list, lambda_center_list
    def doppler(self,lambda_0,llambda):
        v_r = self.c*(llambda-lambda_0)/lambda_0
        return v_r
    def plot_vr(self):
        lambda_list = []
        lambda_0 = 656.30   #[nm]
        for i in range(len(self.days)):
            F_min_list, sigma_list, lambda_center_list = self.auto_find(i)
            lambda_list.append(self.least_squares(F_min_list,sigma_list,lambda_center_list,i)[2])
        v_array = self.doppler(lambda_0, np.array(lambda_list))
        v_array2 = v_array-(np.sum(v_array)/len(v_array))
        plt.plot(self.days,v_array2)


D6().plot_vr()
plt.xlabel("tid[dager]")
plt.ylabel("hastighet[m/s]")
plt.title("Hastighet uten egenhastighet")
plt.tight_layout()
plt.legend()
plt.show()

for i in range(10):
    F_min_list, sigma_list, lambda_center_list = D6().auto_find(i)
    F_min, sigma, lambda_center = D6().least_squares(F_min_list, sigma_list, lambda_center_list, i)
    D6().plot_data(i)
    F_model = D6().F_model(1,F_min, sigma, lambda_center, D6().get_data()[i])
    plt.plot(D6().get_data()[1,:,0],F_model,label="Least squares")
    plt.legend()
    plt.show()
F_min, sigma, lambda_center = D6().least_squares(F_min_list = [0.8,0.7], sigma_list = [0.001,0.02], lambda_center_list = [656.35,656.37], x = 1)
D6().plot_data(1)
plt.title("Minste kvadraters med øyemål")
F_model = D6().F_model(1, F_min, sigma, lambda_center, D6().get_data()[1])
plt.plot(D6().get_data()[1,:,0],F_model)      

plt.show()

for i in range(10):
    D6().plot_data(x = i)
    plt.show()

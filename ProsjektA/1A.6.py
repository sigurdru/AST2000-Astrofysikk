import numpy as np
import random as random
import ast2000tools.utils as utils
seed = utils.get_seed('sigurdru')
np.random.seed(seed)
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as constants
system = SolarSystem(seed)

class rocket():
    def __init__(self,n):
        self.solmasse = 1.98892*10**30                                  #[kg] solmasse i følge wikipedia
        self.G = constants.G                                            #[m^3/kg/s^2] gravitasjonskonstanten til planeten min
        self.n = n                                                      #antall partikler
        self.L = 10**-6                                                 #[m] lengden til en av sidekantene til boksen
        self.k = constants.k_B                                          #[m^2*kg/s^2/K] Boltmanns konstant
        self.T = 10000                                                  #[K] temperaturen til gassen
        self.m_H2 = constants.m_H2                                      #[kg] massen til H2 molekylet
        self.sigma = np.sqrt(self.k*self.T/self.m_H2)                   #standard avvik
        self.m_satellite = 1000                                         #[kg] massen til satelitten
    def positions(self):
        #finner en nx3 matrise med tilfeldige posisjoner.
        positions = np.random.uniform(low=-self.L/2, high=self.L/2, size=(self.n,3))
        return positions
    def velocities(self):
        #finner en nx3 matrise med hastigheter som følger normalfordelingen.
        velocities = np.random.normal(0,self.sigma,size=(self.n,3))
        return velocities
    def mean_ke(self):
        #regner ut forventet og faktisk gjennomsnittlig kinetisk energi.
        expected_K = (3/2)*self.k*self.T
        vel = rocket(self.n).velocities()
        actual_K = np.sum(np.linalg.norm(vel,axis = 1)**2)*0.5*self.m_H2/self.n
        return expected_K,actual_K
    def mean_vel(self):
        #regner ut forventet og faktisk gjenomsnittlig fart.
        expected_mean_vel = 2*np.sqrt((2*self.k*self.T)/(np.pi*self.m_H2))
        vel = rocket(self.n).velocities()
        actual_mean_vel = np.sum(np.linalg.norm(vel,axis=1))/self.n
        return expected_mean_vel, actual_mean_vel
    def simulation(self):
        #returnerer forventet trykk, faktisk trykk, antall partikler som rømmer, endring i bevegelsesmengde
        dt = 10**-12                            #tidssteg
        steps = 1000                            #antall tidssteg
        areal = self.L**2*6                     #areal til boksen
        vel = rocket(self.n).velocities()       #henter hastighet- og posisjons-array
        pos = rocket(self.n).positions()
        particles_escaping = 0                  #lager en variabel som teller antall partikler som rømmer
        speed = 0                               #variabel som teller hastighet fluksen ut av boksen
        for _ in range(steps):                  #løper gjennom tidsstegene, oppdaterer posisjon og registrerer kollisjoner
            pos += vel*dt
            index = np.where(abs(pos)>self.L/2) #finner indekseringen til partikler utafor boksen.
            speed += np.sum(abs(vel[index]))
            vel[index] *= (-1)                  #snur hastighetskomponenten til partiklene som kolliderer.
            particles_escaping += len(index[0]) #teller hvor mange partikler som kolliderer med veggen.
        fuel_loss = particles_escaping*self.m_H2/24 #hvor mye drivstoff som blir brukt
        expected_pressure = (self.n/self.L**3)*self.k*self.T    #forventet trykk
        force = 2*speed*self.m_H2/(steps*dt)    #kraften som partiklene gjør på veggen
        change_in_momentum = speed*self.m_H2/24 #bevegelsesmengde fluks ut av hullet
        actual_pressure = force/areal           #regner ut faktisk trykk
        return expected_pressure, actual_pressure, fuel_loss, change_in_momentum
    
    def acceleration(self):
        #regner ut hvor mange rakettbokser vi trenger og hvor mye drivstoff
        fuel_loss, change_in_momentum = rocket(self.n).simulation()[2:] #henter hvor mye drivstoff en boks bruker, og endring i bevegelsesmengde
        akselerasjon = change_in_momentum*10**9/self.m_satellite    #regner ut hvor mye en boks akselererer raketten
        masse_0 = system.masses[0]*self.solmasse                        #[kg] massen til "home planet"
        radius_0 = 1000*system.radii[0]                        #[m] radius til "home planet" (ganger med 1000 siden enhet er originalt km)
        time = 20*60                #tiden vi ønsker å bruke
        boxes = np.roots([akselerasjon**3*0.5*time**4,akselerasjon**2*radius_0*time**2,0,-2*self.G*masse_0])[-1] #finner hvor mange bokser vi trenger
        fuel_loss = fuel_loss*boxes*time*10**9 #hvor mye drivstoff vi bruker
        v_esc = akselerasjon*boxes*time #hva rømningshastigheten er.

        print(v_esc)
        print(boxes)
        print(fuel_loss)


        

expected_K, actual_K = rocket(n = 10**5).mean_ke()
print("|ForventetE_k: %e|FaktiskE_k: %e|forhold: %f" \
               %(expected_K,actual_K,expected_K/actual_K))

expected_mean_vel, actual_mean_vel = rocket(n = 10**5).mean_vel()
print("|Forventet fart: %f|faktisk fart: %f|forhold: %f" \
              %(expected_mean_vel,actual_mean_vel,expected_mean_vel/actual_mean_vel))

expected_pressure, actual_pressure =rocket(n = 10**5).simulation()[0:2]
print("|forventet trykk: %f|faktsik trykk: %f|forhold: %f" \
        %(expected_pressure,actual_pressure,expected_pressure/actual_pressure))

rocket(n=10**5).acceleration()

import numpy as np
import matplotlib.pyplot as plt
import sys

def P(T):
    mean = 29
    st_dev = 2
    P = (1/(st_dev*np.sqrt(2*np.pi)))*np.exp(-(((T-mean)**2)/(2*st_dev**2)))
    #P = 1/(np.exp((T-mean)**2/2*st_dev**2)*(st_dev*np.sqrt(2*np.pi)))
    return P

x_array = np.linspace(20,40,1000)
y_array = P(x_array)
plt.plot(x_array,y_array)
#plt.show()

a = 27
b = 31
n = 1000
dT = (b-a)/n
#oppgave 1 trapesmetoden
temp = np.linspace(a,b,n)
"""prob = dT/2*(P(temp[1]) + P(temp[-1]))
for i in temp[1:-1]:
    prob += dT*P(i)"""
prob = 0
for i in range(n):
    prob += P(a+i*dT)*dT

print(prob)
#oppgave 2
"""
sannsynligheten er 68% siden temperaturen ligger mellom 27 og 31 som er mean +- standard deviation
"""
#oppgave 3
#bruker resultatet fra oppgave 1 og ganger det med seg selv 7 ganger, da får jeg resultatet
"""
for i in range(7):
    prob *= prob
print (prob)"""
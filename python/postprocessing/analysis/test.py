from __future__ import division
import matplotlib.pyplot as plt
import numpy as np


def AsimovSignificance(s,b):
    #print(2*(s+b)*np.log(s/b)-s)
    return np.sqrt(2*((s+b)*np.log(1+s/b)-s))

#za2 = 2.*( (s+b) * std::log(1. + s/b) -s );
#sqrt(za2)

# signal1 = np.array([30, 100, 200, 500, 600, 800, 1000, 900, 600, 500, 300, 200, 10])
# signal2 = signal1*10
# bckgrd1 = 10000*np.array([50, 55, 54, 40, 49, 59, 45, 56, 30, 50 ,40 ,10 ,20])
# bckgrd2 = bckgrd1/20 +  10*np.random.rand(13)

graph1 = []
graph2 = []

ratio1 = [0.1, 250]
ratio2 = [0.01, 2]


x = np.linspace(0, 500, 1000)

for i in x:
    graph1.append(AsimovSignificance(i*ratio1[0], ratio1[1]))
    graph2.append(AsimovSignificance(i*ratio2[0], ratio2[1]))
    
    
    
    # final1 = 0
    # final2 = 0
    # for j in range(len(signal1)):
    #     final1+=AsimovSignificance(i*signal1[j], bckgrd1[j])**2
    #     final2+=AsimovSignificance(i*signal2[j], bckgrd2[j])**2
    # final1 = np.sqrt(final1)
    # final2 = np.sqrt(final2)
    # graph1.append(final1)
    # graph2.append(final2)

plt.plot(x, graph1)
plt.plot(x,graph2)
plt.show()


# for i in range(len(signal1)):
#     final1
# x = np.linspace(0, 400000, 1000)
# plt.plot(x, AsimovSignificance(x, 5000))
# plt.plot(x, AsimovSignificance(x, 0.1))
# plt.show()
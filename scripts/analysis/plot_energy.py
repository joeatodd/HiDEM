import numpy as np
import matplotlib.pyplot as plt

dtopr = np.genfromtxt("dtopr", delimiter=", ")
dtop00 = np.genfromtxt("dtop00")
dtop01 = np.genfromtxt("dtop01")


fig_max = np.max(dtop00[:,(1,2,4)])
fig_min = np.min(dtop00[:,(1,2,4)])
buff = (fig_max - fig_min) * 0.05
fig_max+=buff
fig_min-=buff

fig, ax1 = plt.subplots()
ax1.plot(dtop00[:,0],dtop00[:,1],label="wens")
ax1.plot(dtop00[:,0],dtop00[:,2], label="elast")
ax1.plot(dtop00[:,0],dtop00[:,4], label="pot")
ax1.set_ylim((fig_min, fig_max))
ax1.plot(dtop00[:,0],dtop00[:,3], label="kins")
plt.legend(loc=0)

fig, ax1 = plt.subplots()
ax1.plot(dtop01[:,0],dtop01[:,1],label="drag")
ax1.plot(dtop01[:,0],dtop01[:,2],label="damp")
ax1.plot(dtop01[:,0],dtop01[:,3],label="pressure")
ax1.plot(dtop01[:,0],dtop01[:,4],label="bed interact")


plt.plot(data[:,0], data[:,1])

plt.savefig("dtopr.png")

import numpy as np
import os, matplotlib, time, glob
import matplotlib.pyplot as plt



# get the data from txt
path = os.path.dirname(os.path.realpath(__file__))
savePath = f"{path}/1200/0.013-0.500-0.500/90.00-44.55"
spinAngles = np.genfromtxt(f"{savePath}/spinAngles.txt")
posArray = np.genfromtxt(f"{path}/1200/system.txt")
maxZ = np.where(posArray[:,2] == np.amax(posArray[:,2]))[0]
posArray = posArray[maxZ,0:2]
hyst = np.genfromtxt(f"{savePath}/hysteresis_4.txt")
# files = glob.glob(f"{savePath}/*Point*.txt")
# points = []
# for el in files:
#     file = el.split("/")[-1].split("_")
#     name = file[0]
#     value = file[1].split(".txt")[0]
#     points.append((name,value))

fig, ax = plt.subplots(figsize=(12, 10))
# fig, (ax1,ax2) = plt.subplots(1,2,figsize=(22, 10),gridspec_kw={'width_ratios': [1.26, 1]}) # corrige os tamanhos
# plt.tight_layout()
# ax1.set_xlim((0,20))
# ax2.set_ylim((-1,1))

norm = matplotlib.colors.Normalize(vmin=-1,vmax=1)
cm = matplotlib.cm.turbo
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])
cm2 = matplotlib.cm.seismic
# plt.colorbar(sm,location="bottom")
fig.colorbar(sm,ax=ax)
ax.set_box_aspect(1)
magneticFieldStep = 0.005 # testando passo
maxExternalMagneticField = 4
numberOfSteps = int(maxExternalMagneticField/magneticFieldStep)
magneticField = [maxExternalMagneticField]

X = np.arange(0,20.5,1)
Y = np.arange(0,20.5,1)
x,y = np.meshgrid(X,Y)

h = [4.000,0.300,0.000,-0.260,-0.290,-0.320,-0.450]
# h = [4.000,1.000,0.000,-0.270,-0.320,-0.340,-0.500,-0.600,-0.750,-1.520,0.350]
hIndexes = []
# add the values to the list, completing the loop to the hysteresis curve
for i in range(1,4*numberOfSteps+1):
    step = -magneticFieldStep if 2*numberOfSteps + 1 > i >= 1 else magneticFieldStep
    magneticField.append(magneticField[-1] + step)
for p in h:
    indexes = np.where(np.round(hyst[:,0],3) == p)[0]
    for l in indexes:
        hIndexes.append(l)
# input(hIndexes)
# plt.tick_params(axis='y', which='both', labelleft='on', labelright='on',labeltop="on")
# plt.tick_params(axis='x', which='both',labeltop="on")
# ax.set_axis_off()
# ax.grid()
# ax.plot(hyst[:,0],hyst[:,1],marker="o",markersize=2,color="red",linewidth="1")
# for point,value in points:
for value in hIndexes:
    ax.clear()
    posDataArray = spinAngles[value,:].reshape(-1,3)
    posDataArray = np.c_[posArray[:,0],posArray[:,1],posDataArray]
    theta = posDataArray[:,2]
    phi = posDataArray[:,3]
    mag = round(hyst[value,1],2)
    h = round(hyst[value,0],2)
    ax.set_title(f"h = {h}")
    # Z = posDataArray[:,4].reshape(20,20)
    Z = np.cos(theta).reshape(20,20)
    Z2 = np.c_[Z[:,0],Z,Z[:,-1]]
    Z2 = np.r_[[Z2[0,:]],Z2,[Z2[-1,:]]]
    d = ax.contourf(Z2,cmap="turbo",norm=norm,levels=100,origin=None,extent=[0,20,0,20])
    ax.quiver(posDataArray[:,0],posDataArray[:,1],np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),scale=35,color="white",headwidth=6,headlength=6)
    for c in d.collections:
        c.set_edgecolor("face")
    plt.savefig(f"points45/{value}.pdf",bbox_inches="tight")
# plt.savefig(f"{savePath}/hysteresis.pdf",bbox_inches="tight")
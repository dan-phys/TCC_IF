import numpy as np
import os, matplotlib, time
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from matplotlib.backends.backend_pdf import PdfPages

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Evolução do sistema magnético', artist='Daniel',
                comment='Usado para o TCC')
fps = 30
writer = FFMpegWriter(fps=fps, metadata=metadata)

# get the data from txt
path = os.path.dirname(os.path.realpath(__file__))
savePath = f"{path}/1200/0.013-0.500-0.500/90.00-44.55"
posArray = np.genfromtxt(f"{path}/1200/system.txt")
hyst = np.genfromtxt(f"{savePath}/hysteresis_4.txt")
maxZ = np.where(posArray[:,2] == np.amax(posArray[:,2]))[0]
posDataArray = posArray[maxZ,0:2]
dataArray = np.genfromtxt(f"{savePath}/spinAngles.txt")
# maxColumns = dataArray.shape[1]
# allTheta = np.genfromtxt(f"{savePath}/spinAngles.txt",usecols=([k for k in range(0,maxColumns,3)]))

# fig, (ax1,ax2) = plt.subplots(1,2,figsize=(15, 8))
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(22, 10),gridspec_kw={'width_ratios': [1.26, 1]}) # corrige os tamanhos
plt.tight_layout()
# ax1.set_xlim((0,20))
ax2.set_ylim((-1,1))
norm = matplotlib.colors.Normalize(vmin=-1,vmax=1)
# norm.autoscale(np.cos(allTheta))
cm = matplotlib.cm.turbo
cm2 = matplotlib.cm.seismic
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])
# plt.colorbar(sm,location="bottom")
fig.colorbar(sm,ax=ax1)
ax1.set_box_aspect(1)
ax2.set_box_aspect(1)
start = time.time()
magneticFieldStep = 0.005 # testando passo
maxExternalMagneticField = 4
numberOfSteps = int(maxExternalMagneticField/magneticFieldStep)
magneticField = [maxExternalMagneticField]
X = np.arange(0,20.5,1)
Y = np.arange(0,20.5,1)
x,y = np.meshgrid(X,Y)

# add the values to the list, completing the loop to the hysteresis curve
for i in range(1,4*numberOfSteps+1):
    step = -magneticFieldStep if 2*numberOfSteps + 1 > i >= 1 else magneticFieldStep
    magneticField.append(magneticField[-1] + step)

h = 0
# hyst[:,0] *= (4*np.pi*1710)
with writer.saving(fig, f"{savePath}/spinAngles.mp4", 200):
    for line in dataArray:
        ax1.clear()
        ax2.clear()
        array = line.reshape(-1,3)
        theta = array[:,0]
        phi = array[:,1]
        magField = hyst[0:h,0]
        mag = hyst[0:h,1]
        ax1.set_title(f"h = {magneticField[h]:.2f}")
        Z = array[:,2].reshape(20,20)
        Z2 = np.c_[Z[:,0],Z,Z[:,-1]]
        Z2 = np.r_[[Z2[0,:]],Z2,[Z2[-1,:]]]
        ax1.contourf(Z2,cmap="turbo",norm=norm,levels=100,origin=None,extent=[0,20,0,20])
        ax1.quiver(posDataArray[:,0],posDataArray[:,1],np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),scale=35,color="white",headwidth=6,headlength=6)
        ax2.grid()
        ax2.plot(hyst[0,0],hyst[0,1],marker="o",markersize=0.1,alpha=0)
        ax2.plot(hyst[1600,0],hyst[1600,1],marker="o",markersize=0.1,alpha=0)
        ax2.set_title("Histerese")
        ax2.plot(magField,mag,marker="o",markersize=2,color="red",linewidth="1")
        writer.grab_frame()
        h += 1
        progress = h/len(magneticField)*100
        print(f"Progresso: [----- {progress:.2f}% -----]", end= "\r" if progress < 100 else "\n")
end = time.time()
print(end - start)
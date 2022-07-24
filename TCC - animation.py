import numpy as np
import os, matplotlib, time
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Evolução do sistema magnético', artist='Daniel',
                comment='Usado para o TCC')
fps = 30
writer = FFMpegWriter(fps=fps, metadata=metadata)

# get the data from txt
path = os.path.dirname(os.path.realpath(__file__))
savePath = f"{path}/1200/0.013-0.500-0.500/90.00-45.00"
posArray = np.genfromtxt(f"{path}/1200/system.txt")
maxZ = np.where(posArray[:,2] == np.amax(posArray[:,2]))[0]
posDataArray = posArray[maxZ,0:2]
dataArray = np.genfromtxt(f"{savePath}/spinAngles.txt")


fig, ax = plt.subplots()
plt.axis('scaled')
plt.xlim((0,20))
plt.ylim((0,20.5))
ax.set_yticks([])
ax.set_xticks([])
norm = matplotlib.colors.Normalize(vmin=-1,vmax=1)
cm = matplotlib.cm.turbo
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])
plt.colorbar(sm)
start = time.time()
magneticFieldStep = 0.005 # testando passo
maxExternalMagneticField = 4
numberOfSteps = int(maxExternalMagneticField/magneticFieldStep)
magneticField = [maxExternalMagneticField]

# add the values to the list, completing the loop to the hysteresis curve
for i in range(1,4*numberOfSteps+1):
    step = -magneticFieldStep if 2*numberOfSteps + 1 > i >= 1 else magneticFieldStep
    magneticField.append(magneticField[-1] + step)

h = 0
with writer.saving(fig, f"{savePath}/spinAngles.mp4", 200):
    for line in dataArray:
        ax.clear()
        array = line.reshape(-1,3)
        theta = array[:,0]
        phi = array[:,1]
        ax.set_title(f"h = {magneticField[h]:.2f}")
        ax.quiver(posDataArray[:,0],posDataArray[:,1],np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),scale=25,color=cm(norm(array[:,2])))
        writer.grab_frame()
        h += 1
        progress = h/len(magneticField)*100
        print(f"Progresso: [----- {progress:.2f}% -----]", end= "\r" if progress < 100 else "\n")
end = time.time()
print(end - start)
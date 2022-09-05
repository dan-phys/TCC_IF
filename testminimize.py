from unittest import result
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Diferença entre os métodos de minimização', artist='Daniel',
                comment='Usado para o TCC')
fps = 1
writer = FFMpegWriter(fps=fps, metadata=metadata)

count = 0

bfgs = []
stepest_descend = []

def energy(param):
    t,p = param
    mz = np.cos(t)
    my = np.sin(t)*np.sin(p)
    mx = np.sin(t)*np.cos(p)
    aniso = (mx**2*my**2 + mx**2*mz**2 + my**2*mz**2)
    return aniso

def callback(xk):
    global count
    bfgs.extend([xk[0],xk[1]])
    count += 1

angles = (80*np.pi/180,60*np.pi/180)
bfgs.extend(angles)
stepest_descend.extend(angles)
resultMin = minimize(energy,angles,tol=1e-4,method="BFGS",callback=callback)

a = 0.1
n = 0
x,y = angles
mx = np.sin(x)*np.cos(y)
my = np.sin(x)*np.sin(y)
mz = np.cos(x)        
dmxt = np.cos(x)*np.cos(y)
dmxf = -np.sin(x)*np.sin(y)
dmyt = np.cos(x)*np.sin(y)
dmyf = np.sin(x)*np.cos(y)
dmzt = -np.sin(x)
dmzf = 0.00000000000
deut=-(mx**3*dmxt + my**3*dmyt + mz**3*dmzt)
deuf=-(mx**3*dmxf + my**3*dmyf + mz**3*dmzf)
gt =  deut
gf = deuf
g = np.sqrt((gt)**2+(gf)**2)
while g >= 10**(-4):
    x = x - gt*a
    y = y - gf*a
    mx = np.sin(x)*np.cos(y)
    my = np.sin(x)*np.sin(y)
    mz = np.cos(x)        
    dmxt = np.cos(x)*np.cos(y)
    dmxf = -np.sin(x)*np.sin(y)
    dmyt = np.cos(x)*np.sin(y)
    dmyf = np.sin(x)*np.cos(y)
    dmzt = -np.sin(x)
    dmzf = 0.00000000000
    deut=-(mx**3*dmxt + my**3*dmyt + mz**3*dmzt)
    deuf=-(mx**3*dmxf + my**3*dmyf + mz**3*dmzf)
    gt =  deut
    gf = deuf
    g = np.sqrt((gt)**2+(gf)**2)
    stepest_descend.extend([x,y])
    n += 1
bfgs = np.array(bfgs).reshape(-1,2)
stepest_descend = np.array(stepest_descend).reshape(-1,2)
# np.savetxt("bfgs.txt",bfgs,fmt='%1.5f',newline="},\n",delimiter="/")
# np.savetxt("stepest_descend.txt",stepest_descend,fmt='%1.5f',newline="},\n",delimiter="/")
# print(bfgs)
# print(stepest_descend)

fig = plt.figure()

ax = fig.add_subplot(111,projection="3d") 
r = [-0.5,0.5]
maxL = 5
X, Y = np.meshgrid(r, r)
one = np.ones(4).reshape(2, 2) - 0.5
x = 0.5
y = 0.5
z = 0.5


c = 1
for line in bfgs:
    ax.clear()
    ax.plot_surface(X + x,Y + y,one + z, color="blue",edgecolors="black",linewidth=0.3,alpha=0.6)
    ax.plot_surface(X + x,Y + y,-one + z, color="blue",edgecolors="black",linewidth=0.3,alpha=0.6)
    ax.plot_surface(X + x,-one + y,Y + z, color="blue",edgecolors="black",linewidth=0.3,alpha=0.6)
    ax.plot_surface(X + x,one + y,Y + z, color="blue",edgecolors="black",linewidth=0.3,alpha=0.6)
    ax.plot_surface(one + x,X + y,Y + z, color="blue",edgecolors="black",linewidth=0.3,alpha=0.6)
    ax.plot_surface(-one + x,X + y,Y + z, color="blue",edgecolors="black",linewidth=0.3,alpha=0.6)
    ax.set_axis_off()
    ax.view_init(azim=-30, elev=30)
    ax.set_title(f"Nº de passos: {c}")
    ax.quiver(x,y,z,np.sin(line[0])*np.cos(line[1]),np.sin(line[0])*np.sin(line[1]),np.cos(line[0]),length = 0.7,color="black")
    plt.savefig(f"bfgs/{c}.png",dpi=200,bbox_inches="tight")
    c += 1
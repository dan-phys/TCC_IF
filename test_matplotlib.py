import matplotlib.pyplot as plt
import numpy as np

partDivision = (20,20,3)
step = 1
maxL = max(partDivision)
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d") 
ax.set_xlim3d(0, maxL)
ax.set_ylim3d(0, maxL)
ax.set_zlim3d(0, maxL)
r = [-0.5,0.5]
X, Y = np.meshgrid(r, r)
one = np.ones(4).reshape(2, 2) - 0.5

# generate the particles positions randomly, according to the model
qntParticles = 0
for k in range(partDivision[2]):
    for j in range(partDivision[1]):
        for i in range(partDivision[0]):
            qntParticles += 1
            x,y,z = (step*np.array([i,j,k]) + 0.5)

            ax.plot_surface(X + x,Y + y,one + z, color="blue",edgecolors="black",linewidth=0.3)
            ax.plot_surface(X + x,Y + y,-one + z, color="blue",edgecolors="black",linewidth=0.3)
            ax.plot_surface(X + x,-one + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)
            ax.plot_surface(X + x,one + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)
            ax.plot_surface(one + x,X + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)
            ax.plot_surface(-one + x,X + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)

ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.zaxis.set_ticklabels([])
color_tuple = (1.0, 1.0, 1.0, 0.0)
ax.w_xaxis.line.set_color(color_tuple)
ax.w_yaxis.line.set_color(color_tuple)
ax.w_zaxis.line.set_color(color_tuple)
ax.tick_params(color=color_tuple, labelcolor=color_tuple)
plt.show()
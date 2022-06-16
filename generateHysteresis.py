import matplotlib.pyplot as plt
import os

path = os.path.dirname(os.path.realpath(__file__))

fig, ax = plt.subplots(figsize=(10,10))
ax.grid(which="both",color="grey", linewidth=1, linestyle="-", alpha=0.6)
ax.set_ylim((-1,1))
ax.set_xlim((-3, 3))
ax.set_xlabel("h", size=14)
ax.set_ylabel("m", size=14)
ax.set_title(f"Curva de histerese: {model} {numOfParticles} part {concentration} conc {intExchange} exc")
plt.plot(magneticField,magnetizationMatrix,color="red",linewidth="1")
fig.savefig(f"{path}/{model}/{numOfParticles}/{concentration}/hysteresis_{intExchange:.1f}_{cutoff}.{ImgFile}",bbox_inches="tight")
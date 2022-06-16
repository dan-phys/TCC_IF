import glob
import matplotlib.pyplot as plt

# fazer para remanência e coercividade

path = "/home/dan_phys/TCC/SC/Hysteresis/Txt"
fig,ax = plt.subplots()
remainingCutoff = []
start = len(path) + 18
end = 9
particles = 13**3
files = glob.glob(f"{path}/hys_{particles}*")
files.sort()
for file in files:
    arq = open(file,"r")
    remaining = arq.readlines()[120].split()
    remainingCutoff.append((float(file[start:-end]),remaining[1]))
    arq.close()

remainingCutoff.sort()
cutoff = []
remaining = []

for data in remainingCutoff:
    cutoff.append(data[0])
    remaining.append(float(data[1]))
print(cutoff)
print(remaining)
plt.plot(cutoff,remaining,color="red",linewidth="1")
fig.savefig(f"testeCutoff_{particles}.png",bbox_inches="tight")

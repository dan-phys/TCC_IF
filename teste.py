import numpy as np
import os
import matplotlib.pyplot as plt

path = os.path.dirname(os.path.realpath(__file__))
savePath = f"{path}/1200/0.013-0.500-0.500/90.00-90.00"
matrix = np.genfromtxt(f"{savePath}/data.txt")

def nextPosition(value):
    if value > 0.10:
        return value/abs(value)*np.ceil(abs(value))
    else:
        return 0

def checkBoundaries(tup:tuple):
    for el in tup:
        if 0 < el < 20:
            continue
        else:
            return
    return tup

def removeNone(lista:list):
    return list(filter(lambda x: type(x) is tuple,lista))

def removeSame(x,y,list):
    pos = (x,y)
    a = []
    for el in list:
        if el == pos:
            continue
        else:
            a.append(el)
    return a


def findVortex(x,y,phi,lista):
    nextX = nextPosition(np.cos(phi))
    nextY = nextPosition(np.sin(phi))
    possiblePositions = [(x + nextX, y), (x + nextX, y + nextY), (x, y + nextY)]
    possiblePositions = removeNone(map(checkBoundaries,possiblePositions))
    possiblePositions = removeSame(x,y,possiblePositions)
    a = [[],[],[]]
    if len(possiblePositions) > 0:
        c = 0
        for pos in possiblePositions:
            if pos in lista[c]:
                break
            lista[c].append(pos)
            partIndex = np.where((matrix[:,0] == pos[0]) & (matrix[:,1] == pos[1]))[0][0]
            b = findVortex(matrix[partIndex,0],matrix[partIndex,1],matrix[partIndex,3],lista)
            c += 1
    #         a[0].append(b[0])
    #         a[1].append(b[1])
    #         a[2].append(b[2])
    # else:
    #     return a

a = [[],[],[]]
x = matrix[0,0]
y = matrix[0,1]
phi = matrix[0,3]
# nextX = nextPosition(np.cos(phi))
# nextY = nextPosition(np.sin(phi))
# possiblePositions = {(x + nextX, y), (x + nextX, y + nextY), (x, y + nextY)}
# possiblePositions = removeNone(map(checkBoundaries,possiblePositions))
# possiblePositions = removeSame(x,y,possiblePositions)
b = findVortex(x,y,phi,a)
fig,(ax1,ax2,ax3) = plt.subplots(1,3)

print(a)
# for lista in a:
# sla = a[0]


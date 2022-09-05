import numpy as np
import matplotlib.pyplot as plt

def findVortex(x,y,phi):
    nextVec = []
    nextX = np.cos(phi)
    nextY = np.sin(phi)
    # olhar para quem tem a maior intensidade, se não olhar para a diagonal
    if abs(nextX) < 0.2:
        X = -1
    if abs(nextY) < 0.2:
        Y = -1
    
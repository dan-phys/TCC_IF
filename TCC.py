import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import os, time

class MagneticSystem:

    def __init__(self, numOfParticles, modelType, concentration,ImgFile):

        # all particles have the same radius
        radius = 1

        # save the data of the system
        self.model = modelType
        self.num = numOfParticles
        self.conc = concentration

        # catch the path to save the image/txt output
        self.path = os.path.dirname(os.path.realpath(__file__))

        # save the positions to verify if the new random 
        # position already have one particle or not
        self.posWithParticles = []

        # calculate the step and the quantity of particler per edge
        # according to the modelType
        if modelType == "SC":
            step = 2*radius
            qntEdge = round(np.cbrt(numOfParticles/concentration))
        elif modelType == "FCC":
            step = radius*np.sqrt(2)
            qntEdge = round(np.cbrt(2*numOfParticles/concentration))
        elif modelType == "BCC":
            step = 2*radius/np.sqrt(3)
            qntEdge = round(np.cbrt(4*numOfParticles/concentration))
        else:
            print("Insira um modelo válido (SC, FCC ou BCC).")
            time.sleep(2)
            exit()
        # the qntEdge must be an integer, so to avoid approximation errors we use round function

        edge = qntEdge*step

        # create the txt file that will save the particles positions
        try:
            file = open(f"{self.path}/{modelType}/{numOfParticles}_{concentration}.txt", "w")
        except FileNotFoundError:
            os.mkdir(f"{self.path}/{modelType}")
            file = open(f"{self.path}/{modelType}/{numOfParticles}_{concentration}.txt", "w")

        # if ImgFile is not False, create an image that represents the particles in space
        if ImgFile:
            fig = plt.figure()
            ax = fig.add_subplot(111,projection="3d")
            u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
            X, Y, Z = radius*np.sin(v)*np.cos(u), radius*np.sin(v)*np.sin(u), radius*np.cos(v)

        # generate the particles positions randomly, according to the model
        qntParticles = 0
        while qntParticles < numOfParticles:
            i,j,k = np.random.randint((qntEdge,qntEdge,qntEdge))
            #x,y,z = coordinates = np.array((i,j,k))*step
            x,y,z = coordinates = [i*step,j*step,k*step]
            if coordinates in self.posWithParticles:
                continue
            else:
                if modelType == "SC":
                    pass
                elif modelType == "FCC":
                    if k%2 == 0:
                        if j%2 == 0:
                            if i%2 == 0:
                                pass
                            else: 
                                continue
                        else:
                            if not i%2 == 0:
                                pass
                            else:
                                continue
                    else:
                        if j%2 == 0:
                            if not i%2 == 0:
                                pass
                            else:
                                continue
                        else:
                            if i%2 == 0:
                                pass
                            else:
                                continue
                elif modelType == "BCC":
                    if k%2 == 0:
                        if j%2 == 0:
                            if i%2 == 0:
                                pass
                            else:
                                continue
                        else:
                            continue
                    else:
                        if not j%2 == 0:
                            if not i%2 == 0:
                                pass
                            else:
                                continue
                        else:
                            continue
                qntParticles += 1
                self.posWithParticles.append(coordinates)
                if qntParticles == 1:
                    pos = f"{x} {y} {z}"
                else:
                    pos = f"\n{x} {y} {z}"
                file.write(pos)

                if ImgFile:
                    ax.plot_surface(X + x, Y + y, Z + z,color="blue")

            progress = qntParticles/numOfParticles*100
            print(f"Calculando as posições: [----- {progress:.0f}% -----]", end= "\r" if progress < 100 else "\n")
        file.close()

        particlesVolume = qntParticles*4*np.pi*radius**3/3
        boxVolume = edge**3
        print(f"Concentração da amostra: {particlesVolume/boxVolume*100:.2f} %")

        if ImgFile:
            print("Gerando imagem...")
            plt.savefig(f"{self.path}/{modelType}/system_{numOfParticles}_{concentration}.{ImgFile}",bbox_inches="tight")

    def minimization(self,maxExternalMagField,intExchange,ImgFile):
        
        # save the positions to use in the minimization
        posValues = {
                        "x":     [],
                        "y":     [],
                        "z":     [], 
                        "phi":   [], 
                        "theta": [], 
                        "Phi":   [], # the phi of the easy magnetization axis
                        "Theta": [], # the theta of the easy magnetization axis
                        "hi_x":  [], 
                        "hi_y":  [], 
                        "hi_z":  []
                    }

        # create the txt file for the hysteresis points
        hysTxt = open(f"{self.path}/{self.model}/hys_{self.num}_{maxExternalMagField}_{intExchange}.txt", "w")

        # if ImgFile is not False, create the image
        if ImgFile:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.grid(which="both",color="grey", linewidth=1, linestryle="-", alpha=0.6)
            ax.set_ylim((-1,1))
            ax.set_xlim((-maxExternalMagField, maxExternalMagField))
            ax.set_xlabel("h", size=14)
            ax.set_ylabel("m", size=14)
            ax.set_title(f"Curva de histerese: {self.model} {self.num} part {self.conc} conc")

        


system = MagneticSystem(15**3/2,"FCC",1,"png")
# system.minimization()
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import os, time
from threading import Thread

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

        # save the positions to use in the minimization
        self.posValues = {
                        "x":     [],
                        "y":     [],
                        "z":     [], 
                        "phi":   [], 
                        "theta": [], 
                        "Phi":   [], # the phi of the easy magnetization axis
                        "Theta": [], # the theta of the easy magnetization axis
                        "hi_x":  [0 for k in range(self.num)], 
                        "hi_y":  [0 for k in range(self.num)], 
                        "hi_z":  [0 for k in range(self.num)]
                    }

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
                self.posValues["x"].append(x)
                self.posValues["y"].append(y)
                self.posValues["z"].append(z)
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

    def minimization(self,maxExternalMagField,intExchange):
        
        self.maxExternalMagField = maxExternalMagField
        self.intExchange = intExchange
        
        magneticFieldStep = 0.05
        numberOfSteps = int(maxExternalMagField/magneticFieldStep)
        self.magneticField = [0]
        minimizationRate = 0.1

        self.magnetizationMatrix = []

        for part in range(self.num):
            phi = np.random.rand()*2*np.pi
            theta = np.random.rand()*np.pi
            Phi = np.random.rand()*2*np.pi
            Theta = np.random.rand()*np.pi
            self.posValues["phi"].append(phi)
            self.posValues["theta"].append(theta)
            self.posValues["Phi"].append(Phi)
            self.posValues["Theta"].append(Theta)

        def firstThreadUseInMinimization(initialRange,finalRange):
            for part in range(initialRange,finalRange):
                particle = [self.posValues["x"][part], self.posValues["y"][part], self.posValues["z"][part]]
                theta = self.posValues["theta"][part]
                phi = self.posValues["phi"][part]
                hi_x = 0.00
                hi_y = 0.00
                hi_z = 0.00
                for k in range(3): # acima, abaixo e na mesma altura da partícula escolhida
                    z = particle[2] + step*(k - 1)
                    for j in range(3):
                        y = particle[1] + step*(j - 1)
                        for i in range(3):
                            x = particle[0] + step*(i - 1)
                            test_particle = [x,y,z]
                            if particle != test_particle and test_particle in self.posWithParticles:
                                dist = np.sqrt((particle[0] - x)**2 + (particle[1] - y)**2 + (particle[2] - z)**2)
                                if dist < 2.1*self.radius:
                                    hi_x += intExchange*np.sin(theta)*np.cos(phi)
                                    hi_y += intExchange*np.sin(theta)*np.sin(phi)
                                    hi_z += intExchange*np.cos(theta)
                self.posValues["hi_x"][part] = hi_x
                self.posValues["hi_y"][part] = hi_y
                self.posValues["hi_z"][part] = hi_z

        def secondThreadUseInMinimization(initialRange,finalRange):
            n = 0
            for part in range(initialRange,finalRange):
                tu = self.posValues["Theta"][part]
                fu = self.posValues["Theta"][part]
                theta = self.posValues["theta"][part]
                phi = self.posValues["phi"][part]
                hi_x = self.posValues["hi_x"][part]
                hi_y = self.posValues["hi_y"][part]
                hi_z = self.posValues["hi_z"][part]

                p3 =  np.cos(theta-tu)-2*np.sin(theta)*np.sin(tu)*np.sin((phi-fu)/2)**2
                dp3t = -(np.sin(theta-tu)+2* np.cos(theta)*np.sin(tu)*np.sin((phi-fu)/2)**2)
                dp3f = -(np.sin(theta)*np.sin(tu)*np.sin(phi-fu))
                dzt = (h+hi_z)*np.sin(theta)-np.cos(theta)*(hi_x*np.cos(phi)+hi_y*np.sin(phi))
                dzf = np.sin(theta)*(hi_x*np.sin(phi)-hi_y*np.cos(phi))
                dut = -p3*dp3t
                duf = -p3*dp3f
                dgt = dzt + dut
                dgf = dzf + duf
                dg = np.sqrt((dgt)**2+(dgf)**2)
                while dg > 0.0001:
                    theta = theta - minimizationRate*dgt
                    phi = phi - minimizationRate*dgf
                    p3 =  np.cos(theta-tu)-2*np.sin(theta)*np.sin(tu)*np.sin((phi-fu)/2)**2
                    dp3t = -(np.sin(theta-tu)+2* np.cos(theta)*np.sin(tu)*np.sin((phi-fu)/2)**2)
                    dp3f = -(np.sin(theta)*np.sin(tu)*np.sin(phi-fu))
                    dzt = (h+hi_z)*np.sin(theta)-np.cos(theta)*(hi_x*np.cos(phi)+hi_y*np.sin(phi))
                    dzf = np.sin(theta)*(hi_x*np.sin(phi)-hi_y*np.cos(phi))
                    dut = -p3*dp3t
                    duf = -p3*dp3f
                    dgt = dzt + dut
                    dgf = dzf + duf
                    dg = np.sqrt((dgt)**2+(dgf)**2)
                self.posValues["theta"][part] = theta
                self.posValues["phi"][part] = phi
                n += 1
            perc = n/(len(self.magneticField)*self.num)*100
            print(f"Calculando a contribuição energética das partículas: [----- {perc:.0f}% -----]", end= "\r" if perc < 100 else "\n")

        # add the values to the list, completing the loop to the hysteresis curve
        for i in range(1,5*numberOfSteps+1):
            step = -magneticFieldStep if 3*numberOfSteps + 1 > i >= numberOfSteps + 1 else magneticFieldStep
            self.magneticField.append(self.magneticField[-1] + step)


        for h in self.magneticField:
            initialRange = 0
            for thread in range(1,10):
                finalRange = int(self.num/9 * thread)
                task = Thread(target=firstThreadUseInMinimization,args=(initialRange,finalRange))
                task.start()
                initialRange = finalRange

            initialRange = 0
            for thread in range(1,10):
                finalRange = int(self.num/9 * thread)
                task = Thread(target=secondThreadUseInMinimization,args=(initialRange,finalRange))
                task.start()
                initialRange = finalRange

            s1 = 0.0
            s2 = 0.0
            for part in range(self.num):
                s1 += np.cos(self.posValues["theta"][part])*np.sin(self.posValues["Theta"][part])
                s2 += np.sin(self.posValues["Theta"][part])
            mag = s1/s2
            self.magnetizationMatrix.append(mag)

        

    def saveHysteresisCurve(self,ImgFile):
        
        # create the txt file for the hysteresis points
        hysTxt = open(f"{self.path}/{self.model}/hys_{self.num}_{self.maxExternalMagField}_{self.intExchange}.txt", "w")
        for el in range(len(self.magneticField)):
            if el == 0:
                hysTxt.write(f"{self.magneticField[el]} {self.magnetizationMatrix[el]}")
            else:
                hysTxt.write(f"\n{self.magneticField[el]} {self.magnetizationMatrix[el]}")
        hysTxt.close()

        # if ImgFile is not False, create the image
        if ImgFile:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.grid(which="both",color="grey", linewidth=1, linestryle="-", alpha=0.6)
            ax.set_ylim((-1,1))
            ax.set_xlim((-self.maxExternalMagField, self.maxExternalMagField))
            ax.set_xlabel("h", size=14)
            ax.set_ylabel("m", size=14)
            ax.set_title(f"Curva de histerese: {self.model} {self.num} part {self.conc} conc")
            plt.plot(self.magneticField,self.magnetizationMatrix,color="red",linewidth="1")
            fig.savefig(f"{self.path}/{self.model}/hys_{self.num}_{self.maxExternalMagField}_{self.intExchange:.1f}.{ImgFile}",bbox_inches="tight")


system = MagneticSystem(15**3,"SC",1,False)
system.minimization(3,0.4)
system.saveHysteresisCurve("png")
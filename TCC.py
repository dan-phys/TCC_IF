import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import os, time
from multiprocessing import Pool

def minimization(test):
    test.minimization()

class MagneticSystem:

    def __init__(self, numOfParticles, modelType, concentration,ImgFile):

        # all particles have the same radius
        self.radius = 1

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
            step = 2*self.radius
            qntEdge = round(np.cbrt(numOfParticles/concentration))
        elif modelType == "FCC":
            step = self.radius*np.sqrt(2)
            qntEdge = round(np.cbrt(2*numOfParticles/concentration))
        elif modelType == "BCC":
            step = 2*self.radius/np.sqrt(3)
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
            X, Y, Z = self.radius*np.sin(v)*np.cos(u), self.radius*np.sin(v)*np.sin(u), self.radius*np.cos(v)

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

        particlesVolume = qntParticles*4*np.pi*self.radius**3/3
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

        indexOfParticles = []

        # add the values to the dict posValues
        for part in range(self.num):
            phi = np.random.rand()*2*np.pi
            theta = np.random.rand()*np.pi
            Phi = np.random.rand()*2*np.pi
            Theta = np.random.rand()*np.pi
            self.posValues["phi"].append(phi)
            self.posValues["theta"].append(theta)
            self.posValues["Phi"].append(Phi)
            self.posValues["Theta"].append(Theta)

        # add the values to the list, completing the loop to the hysteresis curve
        for i in range(1,5*numberOfSteps+1):
            step = -magneticFieldStep if 3*numberOfSteps + 1 > i >= numberOfSteps + 1 else magneticFieldStep
            self.magneticField.append(self.magneticField[-1] + step)

        # separate the indexes for the Pool
        initialRange = 0
        for i in range(1,10):
            finalRange = int(self.num/9 * i)
            tempList = [k for k in range(initialRange,finalRange)]
            initialRange = finalRange
            indexOfParticles.append(tempList)

        indexOfParticles = tuple(indexOfParticles)

        def firstThreadUseInMinimization(indexes):
            for part in range(indexes[0],indexes[-1]+1):
                particle = [self.posValues["x"][part], self.posValues["y"][part], self.posValues["z"][part]]
                theta = self.posValues["theta"][part]
                phi = self.posValues["phi"][part]
                hi_x = 0.00
                hi_y = 0.00
                hi_z = 0.00
                for part2 in range(self.num):
                    x_2 = self.posValues["x"][part2]
                    y_2 = self.posValues["y"][part2]
                    z_2 = self.posValues['z'][part2]
                    theta_2 = self.posValues["theta"][part2]
                    phi_2 = self.posValues["phi"][part2]
                    dist_x = particle[0] - x_2
                    dist_y = particle[1] - y_2
                    dist_z = particle[2] - z_2
                    dist = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    mu_x = np.sin(theta_2)*np.cos(phi_2)*dist_x
                    mu_y = np.sin(theta_2)*np.sin(phi_2)*dist_y
                    mu_z = np.cos(theta_2)*dist_z
                    mu = mu_x + mu_y + mu_z
                    dp = 10**(-3)
                    if part2 != part:
                        hi_x += dp*(3*mu*dist_x/dist**5 - np.sin(theta_2)*np.cos(phi_2)/dist**3)
                        hi_y += dp*(3*mu*dist_y/dist**5 - np.sin(theta_2)*np.sin(phi_2)/dist**3)
                        hi_z += dp*(3*mu*dist_z/dist**5 - np.cos(theta_2)/dist**3)
                    if dist < 2.1*self.radius:
                        hi_x += self.intExchange*np.sin(theta)*np.cos(phi)
                        hi_y += self.intExchange*np.sin(theta)*np.sin(phi)
                        hi_z += self.intExchange*np.cos(theta)
                return hi_x,hi_y,hi_z,part
                # self.posValues["hi_x"][part] = hi_x
                # self.posValues["hi_y"][part] = hi_y
                # self.posValues["hi_z"][part] = hi_z

        def secondThreadUseInMinimization(indexes):
            for part in range(indexes[0],indexes[-1]+1):
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
                # self.posValues["theta"][part] = theta
                # self.posValues["phi"][part] = phi
                return theta,phi,part
                # numberOfInteractions[0] += 1
                # perc = numberOfInteractions[0]/(len(self.magneticField)*self.num)*100
                # print(f"Calculando a contribuição energética das partículas: [----- {perc:.0f}% -----]", end= "\r" if perc < 100 else "\n")
                


        for h in self.magneticField:
            
            if __name__ == "__main__":
                # using multiprocessing to work on 10 processes for the firstThreadUseInMinimization
                tasks = Pool(10)
                resultsFirst = tasks.map(firstThreadUseInMinimization,indexOfParticles)
                tasks.close()
                tasks.join()

                self.posValues["hi_x"][resultsFirst[3]] = resultsFirst[0]
                self.posValues["hi_y"][resultsFirst[3]] = resultsFirst[1]
                self.posValues["hi_z"][resultsFirst[3]] = resultsFirst[2]
                
                # using multiprocessing to work on 10 processes for the secondThreadUseInMinimization

                tasks = Pool(10)
                resultsSecond = tasks.map(secondThreadUseInMinimization,indexOfParticles)
                tasks.close()
                tasks.join()

                self.posValues["theta"][resultsSecond[2]] = resultsSecond[0]
                self.posValues["phi"][resultsSecond[2]] = resultsSecond[1]

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


system = MagneticSystem(1000,"SC",1,False)
system.minimization(3,0.4)
system.saveHysteresisCurve("png")
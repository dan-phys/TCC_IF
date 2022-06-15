import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import os, time
from multiprocessing import Pool


class MagneticSystem:

    def __init__(self, savedPositions = True, numOfParticles = 0, modelType = None, concentration = 0, ImgFile = False):

        if not savedPositions:
            # all particles have the same radius, normalized here to 1
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
                exit()
            # the qntEdge must be an integer, so to avoid approximation errors we use round function

            edge = qntEdge*step

            # create the txt file that will save the particles positions
            try:
                file = open(f"{self.path}/{modelType}/{numOfParticles}/{concentration}/system.txt", "w")
            except FileNotFoundError:
                os.makedirs(f"{modelType}/{numOfParticles}/{concentration}")
                file = open(f"{self.path}/{modelType}/{numOfParticles}/{concentration}/system.txt", "w")

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
                    phi = np.random.rand()*2*np.pi
                    theta = np.random.rand()*np.pi
                    Phi = np.random.rand()*2*np.pi
                    Theta = np.random.rand()*np.pi
                    if qntParticles == 1:
                        pos = f"{x} {y} {z} {phi} {theta} {Phi} {Theta}"
                    else:
                        pos = f"\n{x} {y} {z} {phi} {theta} {Phi} {Theta}"
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
                try:
                    plt.savefig(f"{self.path}/{modelType}/{numOfParticles}/{concentration}/system.{ImgFile}",bbox_inches="tight")
                except FileNotFoundError:
                    os.makedirs(f"{modelType}/{numOfParticles}/{concentration}")
                    plt.savefig(f"{self.path}/{modelType}/{numOfParticles}/{concentration}/system.{ImgFile}",bbox_inches="tight")

    @staticmethod
    def hysteresis(savedPositions,cutoff,maxExternalMagField,intExchange,ImgFile):
        global dipoleAndExchangeContributions,minimization
        path = os.path.dirname(os.path.realpath(__file__))
        fileName = savedPositions.split("/")
        model = fileName[0]
        numOfParticles = int(fileName[1])
        concentration = float(fileName[2])
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
        
        with open(f"{path}/{savedPositions}","r") as file:
            for line in file.readlines():
                data = line.split()
                posValues["x"].append(float(data[0]))
                posValues["y"].append(float(data[1]))
                posValues["z"].append(float(data[2]))
                posValues["phi"].append(float(data[3]))
                posValues["theta"].append(float(data[4]))
                posValues["Phi"].append(float(data[5]))
                posValues["Theta"].append(float(data[6]))
                posValues["hi_x"].append(0)
                posValues["hi_y"].append(0)
                posValues["hi_z"].append(0)

        magneticFieldStep = 0.05
        numberOfSteps = int(maxExternalMagField/magneticFieldStep)
        magneticField = [0]
        minimizationRate = 0.1

        magnetizationMatrix = []

        indexOfParticles = []

        closestParticles = []

        # finding the closest particles for each particle using the cut-off
        start = time.time()
        for part in range(numOfParticles):
            tempList = []
            x1 = posValues["x"][part]
            y1 = posValues["y"][part]
            z1 = posValues["z"][part]
            for part2 in range(numOfParticles):
                if part != part2:
                    x2 = posValues["x"][part2]
                    y2 = posValues["y"][part2]
                    z2 = posValues["z"][part2]
                    dist = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
                    if dist < cutoff:
                        tempList.append(part2)
            closestParticles.append(tempList)   

        print(f"Tempo gasto na procura das partículas mais próximas: {time.time() - start :.1f}s")  

        # add the values to the list, completing the loop to the hysteresis curve
        for i in range(1,5*numberOfSteps+1):
            step = -magneticFieldStep if 3*numberOfSteps + 1 > i >= numberOfSteps + 1 else magneticFieldStep
            magneticField.append(magneticField[-1] + step)

        # separate the indexes for the Pool
        initialRange = 0
        for i in range(1,10):
            finalRange = int(numOfParticles/9 * i)
            tempList = [initialRange,finalRange]
            initialRange = finalRange
            indexOfParticles.append(tempList)

        indexOfParticles = tuple(indexOfParticles)

        def dipoleAndExchangeContributions(indexes):
            hiValues = []
            for part in range(indexes[0],indexes[1]):
                particle = [posValues["x"][part], posValues["y"][part], posValues["z"][part]]
                theta = posValues["theta"][part]
                phi = posValues["phi"][part]
                hi_x = 0.00
                hi_y = 0.00
                hi_z = 0.00
                for part2 in closestParticles[part]:
                    x_2 = posValues["x"][part2]
                    y_2 = posValues["y"][part2]
                    z_2 = posValues['z'][part2]
                    theta_2 = posValues["theta"][part2]
                    phi_2 = posValues["phi"][part2]
                    dist_x = particle[0] - x_2
                    dist_y = particle[1] - y_2
                    dist_z = particle[2] - z_2
                    dist = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    mu_x = np.sin(theta_2)*np.cos(phi_2)*dist_x
                    mu_y = np.sin(theta_2)*np.sin(phi_2)*dist_y
                    mu_z = np.cos(theta_2)*dist_z
                    mu = mu_x + mu_y + mu_z
                    dp = 0.31609375
                    hi_x += dp*(3*mu*dist_x/dist**5 - np.sin(theta_2)*np.cos(phi_2)/dist**3)
                    hi_y += dp*(3*mu*dist_y/dist**5 - np.sin(theta_2)*np.sin(phi_2)/dist**3)
                    hi_z += dp*(3*mu*dist_z/dist**5 - np.cos(theta_2)/dist**3)
                    if dist < 2.1:
                        hi_x += intExchange*np.sin(theta)*np.cos(phi)
                        hi_y += intExchange*np.sin(theta)*np.sin(phi)
                        hi_z += intExchange*np.cos(theta)
                hiValues.append((hi_x,hi_y,hi_z,part))
            return hiValues

        def minimization(indexes):
            anglesValues = []
            for part in range(indexes[0],indexes[1]):
                tu = posValues["Theta"][part]
                fu = posValues["Phi"][part]
                theta = posValues["theta"][part]
                phi = posValues["phi"][part]
                hi_x = posValues["hi_x"][part]
                hi_y = posValues["hi_y"][part]
                hi_z = posValues["hi_z"][part]

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
                anglesValues.append((theta,phi,part))
            return anglesValues

        totalTime = 0
        n = 0
        for h in magneticField:
            start = time.time()

            # using multiprocessing with 9 processes to calculate the contributions to the magnetic field
            tasks = Pool(9)
            resultsFirst = tasks.map(dipoleAndExchangeContributions,indexOfParticles)
            tasks.close()
            tasks.join()

            firstFinal = time.time()

            # update the data
            results = []
            for i in resultsFirst:
                results.extend(i)
            for i in results:
                posValues["hi_x"][i[3]] = i[0]
                posValues["hi_y"][i[3]] = i[1]
                posValues["hi_z"][i[3]] = i[2]
            
            secondStart = time.time()

            # using multiprocessing with 9 processes to minimize the energy
            tasks = Pool(9)
            resultsSecond = tasks.map(minimization,indexOfParticles)
            tasks.close()
            tasks.join()

            secondFinal = time.time()

            # updating the data
            results = []
            for i in resultsSecond:
                results.extend(i)
            for i in results:
                posValues["theta"][i[2]] = i[0]
                posValues["phi"][i[2]] = i[1]

            n += 1
            s1 = 0.0
            s2 = 0.0
            for part in range(numOfParticles):
                s1 += np.cos(posValues["theta"][part])*np.sin(posValues["Theta"][part])
                s2 += np.sin(posValues["Theta"][part])
            mag = s1/s2
            magnetizationMatrix.append(mag)
            perc = n/len(magneticField)*100
            final =  time.time()
            totalTime += final - start
            print(f"Progresso: {perc :.2f}%. Tempo total da interação: {final - start :.1f}s. Tempo total: {totalTime :.1f}s. Tempos primeira e segunda pool: {firstFinal - start:.1f}s e {secondFinal - secondStart :.1f}s.", end= "\r" if perc < 100 else "\n")

    
        # create the txt file for the hysteresis points
        hysTxt = open(f"{path}/{model}/{numOfParticles}/{concentration}/hysteresis_{intExchange:.1f}.txt", "w")
        for el in range(len(magneticField)):
            if el == 0:
                hysTxt.write(f"{magneticField[el]} {magnetizationMatrix[el]}")
            else:
                hysTxt.write(f"\n{magneticField[el]} {magnetizationMatrix[el]}")
        hysTxt.close()

        # if ImgFile is not False, create the image
        if ImgFile:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.grid(which="both",color="grey", linewidth=1, linestryle="-", alpha=0.6)
            ax.set_ylim((-1,1))
            ax.set_xlim((-maxExternalMagField, maxExternalMagField))
            ax.set_xlabel("h", size=14)
            ax.set_ylabel("m", size=14)
            ax.set_title(f"Curva de histerese: {model} {numOfParticles} part {concentration} conc {intExchange} exc")
            plt.plot(magneticField,magnetizationMatrix,color="red",linewidth="1")
            fig.savefig(f"{path}/{model}/{numOfParticles}/{concentration}/hysteresis_{intExchange:.1f}.{ImgFile}",bbox_inches="tight")


system = MagneticSystem()
system.hysteresis(f"SC/{8**3}/1/system.txt",10.1,3,0.4,True)

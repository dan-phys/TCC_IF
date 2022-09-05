import numpy as np
import sympy as sp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os, time
from multiprocessing import Pool

# Co e Fe
# Análise da remanência e coercividade 
# Análise da magnetização inicial: derivada é a suscetibilidade inicial
# Produto BH ****

class MagneticSystem:

    def __init__(self, numOfParticles = 0, modelType = None, concentration = 0, ImgFile = False):
        
        # catch the path to save the image/txt output
        self.path = os.path.dirname(os.path.realpath(__file__))

        savedPositions = os.path.exists(f"{self.path}/{modelType}/{numOfParticles}/{concentration}/system.txt")
        
        if not savedPositions:
            # all particles have the same radius, normalized here to 1
            self.radius = 1

            # save the data of the system
            self.model = modelType
            self.num = numOfParticles
            self.conc = concentration

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
                        pos = f"{x} {y} {z} {theta} {phi} {Theta} {Phi}"
                    else:
                        pos = f"\n{x} {y} {z} {theta} {phi} {Theta} {Phi}"
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
    def hysteresis(savedPositions,maxExternalMagField,intExchange,ImgFile):
        global dipoleAndExchangeContributions,minimization
        path = os.path.dirname(os.path.realpath(__file__))
        fileName = savedPositions.split("/")
        model = fileName[0]
        numOfParticles = int(fileName[1])
        concentration = float(fileName[2])
        posValues = np.genfromtxt(savedPositions)
        his = np.zeros((numOfParticles,3))
        posValues = np.concatenate((posValues,his),axis=1)

        dipConst = 0.31525009
        magneticFieldStep = 0.025
        numberOfSteps = int(maxExternalMagField/magneticFieldStep)
        magneticField = [0]

        magnetizationMatrix = []

        indexOfParticles = []

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
                x_1 = posValues[part,0]
                y_1 = posValues[part,1]
                z_1 = posValues[part,2]
                hi_x = 0.00
                hi_y = 0.00
                hi_z = 0.00
                part2 = np.delete(posValues[:,:],part,0)
                x_2 = part2[:,0]
                y_2 = part2[:,1]
                z_2 = part2[:,2]
                dist_x = x_1 - x_2
                dist_y = y_1 - y_2
                dist_z = z_1 - z_2
                dist = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                # tirar dipParticles depois
                dipParticles = np.where(dist < 16.1)
                x_2 = part2[dipParticles,0]
                y_2 = part2[dipParticles,1]
                z_2 = part2[dipParticles,2]
                theta_2 = part2[dipParticles,3]
                phi_2 = part2[dipParticles,4]
                dist_x = x_1 - x_2
                dist_y = y_1 - y_2
                dist_z = z_1 - z_2
                dist = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                mu_x = np.sin(theta_2)*np.cos(phi_2)*dist_x
                mu_y = np.sin(theta_2)*np.sin(phi_2)*dist_y
                mu_z = np.cos(theta_2)*dist_z
                mu = mu_x + mu_y + mu_z
                hi_x = np.sum(dipConst*(3*mu*dist_x/dist**5 - np.sin(theta_2)*np.cos(phi_2)/dist**3))
                hi_y = np.sum(dipConst*(3*mu*dist_y/dist**5 - np.sin(theta_2)*np.sin(phi_2)/dist**3))
                hi_z = np.sum(dipConst*(3*mu*dist_z/dist**5 - np.cos(theta_2)/dist**3))
                exchangeParticles = np.where(dist < 2.01)
                hi_x += intExchange*np.sum(np.sin(theta_2[exchangeParticles])*np.cos(phi_2[exchangeParticles]))
                hi_y += intExchange*np.sum(np.sin(theta_2[exchangeParticles])*np.sin(phi_2[exchangeParticles]))
                hi_z += intExchange*np.sum(np.cos(theta_2[exchangeParticles]))
                hiValues.append((hi_x,hi_y,hi_z,part))
            return hiValues

        def energyFunc(param,*args):
            h,hi_x,hi_y,hi_z,tu,fu = args
            t, p = param
            mz = np.cos(t)
            my = np.sin(t)*np.sin(p)
            mx = np.sin(t)*np.cos(p)
            dipExc = -(mx*(hi_x) + my*(hi_y) + mz*(hi_z + h))
            aniso = -(np.sin(tu)*np.sin(t)*np.cos(p - fu) +np.cos(tu)*np.cos(t))**2/2
            return dipExc + aniso

        def minimization(indexes):
            anglesValues = []
            for part in range(indexes[0],indexes[1]):
                angles = [posValues[part,3], posValues[part,4]]
                tu = posValues[part,5]
                fu = posValues[part,6]
                hi_x = posValues[part,7]
                hi_y = posValues[part,8]
                hi_z = posValues[part,9]
    
                resultMin = minimize(energyFunc,angles,args=(h,hi_x,hi_y,hi_z,tu,fu),tol=1e-4)
                theta,phi = resultMin.x
                anglesValues.append((theta,phi,part))
            return anglesValues

        totalTime = 0
        n = 0
        minimumCoercivity = 0.1
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
                posValues[i[3],7] = i[0]
                posValues[i[3],8] = i[1]
                posValues[i[3],9] = i[2]
            
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
                posValues[i[2],3] = i[0]
                posValues[i[2],4] = i[1]

            n += 1
            s1 = np.sum(np.cos(posValues[:,3])*np.sin(posValues[:,5]))
            s2 = np.sum(np.sin(posValues[:,5]))
            mag = s1/s2
            magnetizationMatrix.append(mag)
            if round(h,3) == 3:
                np.savetxt("spinExamp-sat.txt",posValues)
            if round(h,3) == 0 and n > 1:
                np.savetxt("spinExamp-reman.txt",posValues)
            if abs(mag) < abs(minimumCoercivity) and h < 0:
                minimumCoercivity = mag
                np.savetxt("spinExamp-coerc.txt",posValues)
            if round(h,3) == -3:
                np.savetxt("spinExamp-satInv.txt",posValues)
            perc = n/len(magneticField)*100
            final =  time.time()
            totalTime += final - start
            print(f"Progresso: {perc :.2f}%. Tempo total da interação: {final - start :.1f}s. Tempo total: {totalTime :.1f}s. Tempos primeira e segunda pool: {firstFinal - start:.1f}s e {secondFinal - secondStart :.1f}s.", end= "\r" if perc < 100 else "\n")

    
        # create the txt file for the hysteresis points
        hysTxt = open(f"{path}/{model}/{numOfParticles}/{concentration}/hysteresis_{intExchange:.2f}.txt", "w")
        for el in range(len(magneticField)):
            if el == 0:
                hysTxt.write(f"{magneticField[el]} {magnetizationMatrix[el]}")
            else:
                hysTxt.write(f"\n{magneticField[el]} {magnetizationMatrix[el]}")
        hysTxt.close()

        # if ImgFile is not False, create the image
        if ImgFile:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.grid(which="both",color="grey", linewidth=1, linestyle="-", alpha=0.6)
            ax.set_ylim((-1,1))
            ax.set_xlim((-maxExternalMagField, maxExternalMagField))
            ax.set_xlabel("h", size=14)
            ax.set_ylabel("m", size=14)
            #ax.set_title("Histerese")
            plt.plot(magneticField,magnetizationMatrix,color="red",linewidth="1")
            fig.savefig(f"{path}/{model}/{numOfParticles}/{concentration}/hysteresis_{intExchange:.2f}.{ImgFile}",bbox_inches="tight")

# concentrações 0.01, 0.2, 0.4, 0.6, 0.8 e 1
numOfParticles = 16**3
concentration = 1.0
system = MagneticSystem(numOfParticles,"SC",concentration,"png")
# for i in range(8,9):
system.hysteresis(f"SC/{numOfParticles}/{concentration}/system.txt",3,0.4,"pdf")

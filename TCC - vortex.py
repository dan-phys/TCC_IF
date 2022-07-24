import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os, time
from multiprocessing import Pool

class MagneticSystem:

    def __init__(self, partDivision, ImgFile = False):
        
        self.path = os.path.dirname(os.path.realpath(__file__))

        numOfParticles = partDivision[0]*partDivision[1]*partDivision[2]

        savedPositions = os.path.exists(f"{self.path}/{numOfParticles}/system.txt")
        
        if not savedPositions:
            
            self.radius = 0.5

            step = 2*self.radius

            # create the txt file that will save the particles positions
            try:
                file = open(f"{self.path}/{numOfParticles}/system.txt", "w")
            except FileNotFoundError:
                os.makedirs(f"{numOfParticles}")
                file = open(f"{self.path}/{numOfParticles}/system.txt", "w")

            # if ImgFile is not False, create an image that represents the particles in space
            if ImgFile:
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

                        if qntParticles == 1:
                            pos = f"{x} {y} {z}"
                        else:
                            pos = f"\n{x} {y} {z}"
                        file.write(pos)
                        if ImgFile:
                            ax.plot_surface(X + x,Y + y,one + z, color="blue",edgecolors="black",linewidth=0.3)
                            ax.plot_surface(X + x,Y + y,-one + z, color="blue",edgecolors="black",linewidth=0.3)
                            ax.plot_surface(X + x,-one + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)
                            ax.plot_surface(X + x,one + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)
                            ax.plot_surface(one + x,X + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)
                            ax.plot_surface(-one + x,X + y,Y + z, color="blue",edgecolors="black",linewidth=0.3)

                        progress = qntParticles/numOfParticles*100
                        print(" "*120, end="\r")
                        print(f"Calculando as posições: [----- {progress:.0f}% -----]", end= "\r" if progress < 100 else "\n")
            file.close()

            if ImgFile:
                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.zaxis.set_ticklabels([])
                color_tuple = (1.0, 1.0, 1.0, 0.0)
                ax.w_xaxis.line.set_color(color_tuple)
                ax.w_yaxis.line.set_color(color_tuple)
                ax.w_zaxis.line.set_color(color_tuple)
                ax.tick_params(color=color_tuple, labelcolor=color_tuple)
                print("Gerando imagem...")
                plt.savefig(f"{self.path}/{numOfParticles}/system.{ImgFile}",bbox_inches="tight",dpi=1000)

    @staticmethod
    def hysteresis(savedPositions,maxExternalMagField,th,ph,ImgFile):
        global dipoleAndExchangeContributions,minimization
        path = os.path.dirname(os.path.realpath(__file__))
        numOfParticles = int(savedPositions)
        # dados do coey: pág 242 (transformar para o cgs)
        # expressão do abo2013 para o lex
        Ms = 1710
        intExchange = 0.5
        dipConst = 0.5
        k = 0.013
        thDegree = f"{th*180/np.pi:.2f}"
        phDegree = f"{ph*180/np.pi:.2f}"
        savePath = f"{path}/{savedPositions}/{k:.3f}-{intExchange:.3f}-{dipConst:.3f}/{thDegree}-{phDegree}"
        posValues = np.genfromtxt(f"{path}/{savedPositions}/system.txt")
        posValues = np.concatenate((posValues,(th,ph,1)*np.ones((numOfParticles,3))),axis=1)
        
        hix = np.zeros(numOfParticles)
        hiy = np.zeros(numOfParticles)
        hiz = np.zeros(numOfParticles)
        maxZ = np.amax(posValues[:,2])
        maxZIndex = np.where(posValues[:,2] == maxZ)[0]

        try:
            os.makedirs(savePath)
        except FileExistsError:
            pass
        
        magneticFieldStep = 0.005 # testando passo
        numberOfSteps = int(maxExternalMagField/magneticFieldStep)
        magneticField = [maxExternalMagField]

        magnetizationMatrix = []

        indexOfParticles = []

        # add the values to the list, completing the loop to the hysteresis curve
        for i in range(1,4*numberOfSteps+1):
            step = -magneticFieldStep if 2*numberOfSteps + 1 > i >= 1 else magneticFieldStep
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
                # a linha da partícula 1 não é considerada, visto que uma partícula não produz campo sobre ela mesma
                part2 = np.delete(posValues[:,:],part,0)
                x_2 = part2[:,0]
                y_2 = part2[:,1]
                z_2 = part2[:,2]
                theta_2 = part2[:,3]
                phi_2 = part2[:,4]
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
                exchangeParticles = np.where(dist < 1.01)
                hi_x += intExchange*np.sum(np.sin(theta_2[exchangeParticles])*np.cos(phi_2[exchangeParticles]))
                hi_y += intExchange*np.sum(np.sin(theta_2[exchangeParticles])*np.sin(phi_2[exchangeParticles]))
                hi_z += intExchange*np.sum(np.cos(theta_2[exchangeParticles]))
                hiValues.append((hi_x,hi_y,hi_z,part))
            return hiValues

        def energyFunc(param,*args):
            h,hi_x,hi_y,hi_z,th,ph,k = args
            t, p = param
            mz = np.cos(t)
            my = np.sin(t)*np.sin(p)
            mx = np.sin(t)*np.cos(p)
            h_x = h*np.sin(th)*np.cos(ph)
            h_y = h*np.sin(th)*np.sin(ph)
            h_z = h*np.cos(th)
            dipExc = -(mx*(hi_x + h_x) + my*(hi_y + h_y) + mz*(hi_z + h_z))
            aniso = k*(mx**2*my**2 + mx**2*mz**2 + my**2*mz**2)
            return dipExc + aniso

        def minimization(indexes):
            anglesValues = []
            for part in range(indexes[0],indexes[1]):
                angles = [posValues[part,3], posValues[part,4]]
                hi_x = hix[part]
                hi_y = hiy[part]
                hi_z = hiz[part]
    
                resultMin = minimize(energyFunc,angles,args=(h,hi_x,hi_y,hi_z,th,ph,k),tol=1e-4,method="BFGS")
                theta,phi = resultMin.x
                anglesValues.append((theta,phi,part))
            return anglesValues

        totalTime = 0
        n = 0

        minimumCoercivity = 0.3
        spinAngles = open(f"{savePath}/spinAngles.txt", "w")
        spinAngles.close()
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
                hix[i[3]] = i[0]
                hiy[i[3]] = i[1]
                hiz[i[3]] = i[2]
            
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
                posValues[i[2],5] = np.sin(i[0])*np.cos(i[1])*np.sin(th)*np.cos(ph) + np.sin(i[0])*np.sin(i[1])*np.sin(th)*np.sin(ph) + np.cos(i[0])*np.cos(th)

            s1 = np.sum(np.cos(posValues[:,3]- th)-2*np.sin(posValues[:,3])*np.sin(th)*np.sin((posValues[:,4] - ph)/2)**2)
            mag = s1/numOfParticles
            magnetizationMatrix.append(mag)

            maxZMatrix = posValues[maxZIndex,3:].reshape(1,3*len(maxZIndex))
            with open(f"{savePath}/spinAngles.txt", "ab+") as f:
                np.savetxt(f,maxZMatrix)

            if round(h,6) == 0 and mag > 0:
                np.savetxt(f"{savePath}/remanenceData_{maxExternalMagField}_{h:.4f}.txt",posValues[:,3:])
            elif abs(mag) < abs(minimumCoercivity) and h < 0:
                minimumCoercivity = mag
                np.savetxt(f"{savePath}/coercivityData_{maxExternalMagField}.txt",posValues[:,3:])

            n += 1
            perc = n/len(magneticField)*100
            final =  time.time()
            totalTime += final - start
            print(" "*123,end = "\r")
            print(f"Progresso: {perc :.2f}%. Tempo total da interação: {final - start :.1f}s. Tempo total: {totalTime :.1f}s. Tempos primeira e segunda pool: {firstFinal - start:.1f}s e {secondFinal - secondStart :.1f}s.", end= "\r" if perc < 100 else "\n")
        
        # create the txt file for the hysteresis points
        with open(f"{savePath}/hysteresis_{maxExternalMagField}.txt", "w") as hysTxt:
            for el in range(len(magneticField)):
                if el == 0:
                    hysTxt.write(f"{magneticField[el]} {magnetizationMatrix[el]}")
                else:
                    hysTxt.write(f"\n{magneticField[el]} {magnetizationMatrix[el]}")

        # if ImgFile is not False, create the image
        if ImgFile:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.grid(which="both",color="grey", linewidth=1, linestyle="-", alpha=0.6)
            ax.set_ylim((-1,1))
            ax.set_xlim((-maxExternalMagField, maxExternalMagField))
            ax.set_xlabel("h", size=14)
            ax.set_ylabel("m", size=14)
            ax.set_title(f"Curva de histerese: {numOfParticles} part {intExchange} exc")
            plt.plot(magneticField,magnetizationMatrix,color="red",linewidth="1")
            fig.savefig(f"{savePath}/hysteresis_{maxExternalMagField}.{ImgFile}",bbox_inches="tight")

sysDist = (20,20,3)
numOfParticles = sysDist[0]*sysDist[1]*sysDist[2]
system = MagneticSystem(sysDist,"png")
# for i in range(60,30,-15):
# rodar para calcular o H
system.hysteresis(numOfParticles,4,np.pi/2,(np.pi*45/180)*0.99999,"pdf")
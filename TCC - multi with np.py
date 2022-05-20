import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import os, time
from multiprocessing import Pool

# all particles have the same radius
radius = 1

# system's data
modelType = "SC"
numOfParticles = 900
concentration = 0.8
concentration = 1
modelImgFile = False
minimization = True
if minimization:
    maxExternalMagField = 3
    intExchange = 0.4
    hysImgFile = "pdf"

# catch the path to save the image/txt output
path = os.path.dirname(os.path.realpath(__file__))

# save the positions to use in the minimization.
# the values are: x y z phi theta Phi Theta hi_x hi_y hi_z
posValues = np.zeros((numOfParticles,10))

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


# if ImgFile is not False, create an image that represents the particles in space
if modelImgFile:
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")
    u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
    X, Y, Z = radius*np.sin(v)*np.cos(u), radius*np.sin(v)*np.sin(u), radius*np.cos(v)

# generate the particles positions randomly, according to the model
qntParticles = 0
while qntParticles < numOfParticles:
    i,j,k = np.random.randint((qntEdge,qntEdge,qntEdge))
    x,y,z = coordinates = [i*step,j*step,k*step]
    if coordinates in posValues[:,0:3].tolist():
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
        
        posValues[qntParticles,0] = x
        posValues[qntParticles,1] = y
        posValues[qntParticles,2] = z
        qntParticles += 1
        if modelImgFile:
            ax.plot_surface(X + x, Y + y, Z + z,color="blue")

    progress = qntParticles/numOfParticles*100
    print(f"Calculando as posições: [----- {progress:.0f}% -----]", end= "\r" if progress < 100 else "\n")

if modelImgFile:
    print("Gerando imagem...")
    plt.savefig(f"{path}/{modelType}/system_{numOfParticles}_{concentration}.{modelImgFile}",bbox_inches="tight")

# save the txt file that will save the particles positions
try:
    np.savetxt(f"{path}/{modelType}/{numOfParticles}_{concentration}", posValues)
except FileNotFoundError:
    os.mkdir(f"{path}/{modelType}")
    np.savetxt(f"{path}/{modelType}/{numOfParticles}_{concentration}", posValues)

if minimization:
    if __name__ == "__main__":
        magneticFieldStep = 0.1
        numberOfSteps = int(maxExternalMagField/magneticFieldStep)
        magneticField = [0]
        minimizationRate = 0.1

        magnetizationMatrix = []

        indexOfParticles = []

        # add the values to the list posValues
        for part in range(numOfParticles):
            phi = np.random.rand()*2*np.pi
            theta = np.random.rand()*np.pi
            Phi = np.random.rand()*2*np.pi
            Theta = np.random.rand()*np.pi
            posValues[part,3] = phi
            posValues[part,4] = theta
            posValues[part,5] = Phi
            posValues[part,6] = Theta

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
        # print(indexOfParticles)
        # print(len(posWithParticles))
        # input()
        def firstThreadUseInMinimization(indexes):
            hiValues = []
            for part in range(indexes[0],indexes[1]):
                particle = [posValues[part,0], posValues[part,1], posValues[part,2]]
                theta = posValues[part,3]
                phi = posValues[part,4]
                hi_x = 0.00
                hi_y = 0.00
                hi_z = 0.00
                for part2 in range(numOfParticles):
                    x_2 = posValues[part2,0]
                    y_2 = posValues[part2,1]
                    z_2 = posValues[part2,2]
                    theta_2 = posValues[part2,3]
                    phi_2 = posValues[part2,4]
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
                        if dist < 2.1*radius:
                            hi_x += intExchange*np.sin(theta)*np.cos(phi)
                            hi_y += intExchange*np.sin(theta)*np.sin(phi)
                            hi_z += intExchange*np.cos(theta)
                hiValues.append((hi_x,hi_y,hi_z,part))
            return hiValues
                # posValues["hi_x"][part] = hi_x
                # posValues["hi_y"][part] = hi_y
                # posValues["hi_z"][part] = hi_z

        def secondThreadUseInMinimization(indexes):
            anglesValues = []
            for part in range(indexes[0],indexes[1]):
                tu = posValues[part,6]
                fu = posValues[part,5]
                theta = posValues[part,4]
                phi = posValues[part,3]
                hi_x = posValues[part,7]
                hi_y = posValues[part,8]
                hi_z = posValues[part,9]

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
                while dg > 0.001:
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
            # using multiprocessing to work on 9 processes for the firstThreadUseInMinimization
            tasks = Pool(9)
            resultsFirst = tasks.map(firstThreadUseInMinimization,indexOfParticles)
            tasks.close()
            tasks.join()
            
            # updating the data
            results = []
            for i in resultsFirst:
                results.extend(i)
            for i in results:
                posValues[i[3],7] = i[0]
                posValues[i[3],8] = i[1]
                posValues[i[3],9] = i[2]
            
            # using multiprocessing to work on 9 processes for the secondThreadUseInMinimization

            tasks = Pool(9)
            resultsSecond = tasks.map(secondThreadUseInMinimization,indexOfParticles)
            tasks.close()
            tasks.join()

            # updating the data
            results = []
            for i in resultsSecond:
                results.extend(i)
            for i in results:
                posValues[i[2],4] = i[0]
                posValues[i[2],3] = i[1]

            s1 = 0.0
            s2 = 0.0
            for part in range(numOfParticles):
                s1 += np.cos(posValues[part,4])*np.sin(posValues[part,6])
                s2 += np.sin(posValues[part,6])
            mag = s1/s2
            magnetizationMatrix.append(mag)
            n += 1
            perc = n/len(magneticField)*100
            final =  time.time()
            totalTime += final - start
            print(f"Progresso: {perc :.2f} %. Tempo gasto para essa interação: {final - start :.1f}s. Tempo total de execução: {totalTime :.1f}s", end= "\r" if perc < 100 else "\n")

            perc = n/len(magneticField) * 100
            print(f"Valores atualizados para o campo magnético {h}. Progresso: {perc :.2f}%.", end= "\r" if perc < 100 else "\n")


        # create the txt file for the hysteresis points
        hysTxt = open(f"{path}/{modelType}/hys_{numOfParticles}_{maxExternalMagField}_{intExchange}.txt", "w")
        for el in range(len(magneticField)):
            if el == 0:
                hysTxt.write(f"{magneticField[el]} {magnetizationMatrix[el]}")
            else:
                hysTxt.write(f"\n{magneticField[el]} {magnetizationMatrix[el]}")
        hysTxt.close()

        # if ImgFile is not False, create the image
        if hysImgFile:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.grid(which="both",color="grey", linewidth=1, linestyle="-", alpha=0.6)
            ax.set_ylim((-1,1))
            ax.set_xlim((-maxExternalMagField, maxExternalMagField))
            ax.set_xlabel("h", size=14)
            ax.set_ylabel("m", size=14)
            ax.set_title(f"Curva de histerese: {modelType} {numOfParticles} part {concentration} conc")
            plt.plot(magneticField,magnetizationMatrix,color="red",linewidth="1")
            fig.savefig(f"{path}/{modelType}/hys_{numOfParticles}_{maxExternalMagField}_{intExchange:.1f}.{hysImgFile}",bbox_inches="tight")
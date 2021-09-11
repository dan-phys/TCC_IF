import os
import matplotlib.pyplot as plt
import numpy as np
import time

def cubic_structure(tp,raio,qnt,tp_arq,*args):
    linha = 0
    pos_w_part = []
    path = os.path.dirname(os.path.realpath(__file__))
    arq1 = open(f"{path}/united_particles_{tp}.txt", "w")
    aprox = 0.1
    if tp == "SC":
        step = 2*raio
        qnt_aresta = (qnt)**(1/3) if not args else (qnt/args[0])**(1/3)
        aprox = 0.09
        # v = (qnt_aresta**3)*4*np.pi/((qnt_aresta+))
    elif tp == "FCC":
        step = raio*np.sqrt(2)
        qnt_aresta = (2*qnt)**(1/3) if not args else int((2*qnt/args[0])**(1/3))
    elif tp == "BCC":
        step = 2*raio/np.sqrt(3)
        qnt_aresta = (4*qnt)**(1/3) if not args else int((4*qnt/args[0])**(1/3))
    if not args:
        if qnt_aresta - int(qnt_aresta) <= aprox:
            qnt_aresta = int(qnt_aresta)
        else:
            qnt_aresta = np.ceil(qnt_aresta)
    aresta = qnt_aresta*step
    if tp_arq == "txt":
        arq = open("{}/{}.txt".format(path,tp), "w")
    elif tp_arq != "txt" and tp_arq != "":
        tipo = tp_arq
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
        X, Y, Z = raio*np.cos(u)*np.sin(v), raio*np.sin(u)*np.sin(v), raio*np.cos(v)
        # Mantém a proporcionalidade na imagem final
        ax.plot_surface(0.0001*np.cos(u)*np.sin(v) - step,0.0001*np.sin(u)*np.sin(v) - step,0.0001*np.cos(v) - step, color="red") # xyz=0
        ax.plot_surface(0.0001*np.cos(u)*np.sin(v) + aresta,0.0001*np.sin(u)*np.sin(v),0.0001*np.cos(v), color="red") # x max
        ax.plot_surface(0.0001*np.cos(u)*np.sin(v),0.0001*np.sin(u)*np.sin(v) + aresta,0.0001*np.cos(v), color="red") # y max
        ax.plot_surface(0.0001*np.cos(u)*np.sin(v),0.0001*np.sin(u)*np.sin(v),0.0001*np.cos(v) + aresta, color="red") # z max
    if args:
        prob = args[0]
        if len(args) > 1:
            m_posx = args[1] # matriz posição x
            m_posy = args[2] # matriz posição x
            m_posz = args[3] # matriz posição x
            m_phi = args[4] # matriz para o ângulo phi
            m_theta = args[5] # matriz para o ângulo theta
    part_qnt = 0
    while part_qnt < qnt:
        i,j,k = np.random.randint(qnt_aresta), np.random.randint(qnt_aresta), np.random.randint(qnt_aresta)
        x,y,z = i*step,j*step,k*step
        coordinates = [x,y,z]
        if coordinates in pos_w_part:
            continue
        else:
            if tp == "SC":
                pass
            elif tp == "FCC":
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
                        if not i%2 ==0:
                            pass
                        else:
                            continue
                    else:
                        if i%2 == 0:
                            pass
                        else:
                            continue
            elif tp == "BCC":
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
        if args:
            p = np.random.rand()
            if p <= prob:
                part_qnt += 1
                pos_w_part.append(coordinates)
                if tp_arq == "txt":
                    if linha == 0: 
                        pos = [f"{x:.2f}", " ", f"{y:.2f}", " ", f"{z:.2f}"]
                        linha += 1
                    else:
                        pos = ["\n", f"{x:.2f}", " ", f"{y:.2f}", " ", f"{z:.2f}"]
                    arq.writelines(pos)
                elif tp_arq == "":
                    m_posx.append(x)
                    m_posy.append(y)
                    m_posz.append(z)
                    m_phi.append(np.random.rand()*2*np.pi)
                    m_theta.append(np.random.rand()*np.pi)
                elif tp_arq != "txt" and tp_arq != "":
                    ax.plot_surface(X + x, Y + y, Z + z,color="blue")
        else:
            part_qnt += 1
            pos_w_part.append(coordinates)
            if tp_arq == "txt":
                if linha == 0: 
                    pos = [f"{x:.2f}", " ", f"{y:.2f}", " ", f"{z:.2f}"]
                    linha += 1
                else:
                    pos = ["\n", f"{x:.2f}", " ", f"{y:.2f}", " ", f"{z:.2f}"] 
                arq.writelines(pos)
            else:
                ax.plot_surface(X + x, Y + y, Z + z, color='blue')
        perc = part_qnt/qnt*100
        print(f"Posições: [----- {perc:.0f}% -----]", end= "\r" if perc < 100 else "\n")
    if tp_arq == "txt":
        arq.close()
    elif tp_arq != "txt" and tp_arq != "":
        plt.savefig("{}/{}.{}".format(path,tp,tipo))

    # linha = 0
    # qnt_united = 0
    # for part in pos_w_part:
    #     unt_part = [tuple(part)]
    #     if tp == "SC":
    #         for pos in range(3): # referente às 3 coordenadas
    #             test_particle = [part[0],part[1],part[2]]
    #             coord = part[pos] - step
    #             for c in range(2):
    #                 test_particle[pos] = coord + 2*c*step
    #                 if test_particle in pos_w_part:
    #                     unt_part.append(tuple(test_particle))
    #     elif tp == "FCC":
    #         for k in range(3): # acima, abaixo e na mesma altura da partícula escolhida
    #             test_particle = [part[0],part[1],part[2]]
    #             z = part[2] + step*(k - 1)
    #             for j in range(3):
    #                 y = part[1] + step*(j - 1)

    #     elif tp == "BCC":
    #         k = 0
    #     if len(unt_part) > 1:
    #         if linha == 0:
    #             arq1.write(str(unt_part)[1:-1])
    #             linha += 1
    #         else:
    #             arq1.write("\n" + str(unt_part)[1:-1])
    #         qnt_united += len(unt_part)

    # qnt_united = 0
    # for part1 in range(len(pos_w_part)):
    #     coord1 = pos_w_part[part1]
    #     for part2 in range(part1 + 1, len(pos_w_part)):
    #         coord2 = pos_w_part[part2]
    #         dist = np.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)
    #         if dist < 2.1*raio:
    #             qnt_united += 1
    # end = time.time()
    # arq1.write(str(qnt_united))
    # arq1.close()
    # print(qnt_united)
    v_cubo = aresta**3
    v_part = qnt*4*np.pi*raio**3/3
    print(f"Concentração: {(v_part/v_cubo)*100:.2f}%")
cubic_structure("BCC",1,2000,"png",0.5)
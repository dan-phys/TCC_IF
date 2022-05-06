import os
import matplotlib.pyplot as plt
import numpy as np
import time
import sympy as sp

def cubic_structure(tp,qnt,tp_arq,prob,minimization,h_max,he_x,tp_hys):
    path = os.path.dirname(os.path.realpath(__file__))
    linha = 0
    raio = 1
    pos_w_part = []
    pos_values = {
                    "x": [],
                    "y": [],
                    "z": [], 
                    "phi": [], 
                    "theta": [], 
                    "phi_u": [], 
                    "theta_u": [], 
                    "hi_x": [0 for k in range(qnt)], 
                    "hi_y": [0 for k in range(qnt)], 
                    "hi_z": [0 for k in range(qnt)]
                }
    aprox = 0.09 
    # mag_txt = open(f"mag_{tp}_{qnt}.txt", "w")

    # ===== Define a separação entre as partículas de acordo com o modelo cúbico ===== #
    if tp == "SC":
        step = 2*raio
        qnt_edge = ((2 - prob)*qnt)**(1/3) # corrigir a probabilidade
    elif tp == "FCC":
        step = raio*np.sqrt(2)
        qnt_edge = (2*(2 - prob)*qnt)**(1/3)
    elif tp == "BCC":
        step = 2*raio/np.sqrt(3)
        qnt_edge = (4*(2 - prob)*qnt)**(1/3)
    if qnt_edge - int(qnt_edge) <= aprox:
        qnt_edge = int(qnt_edge)
    else:
        qnt_edge = np.ceil(qnt_edge)
    # ====================================== ## ======================================= #        
    
    edge = qnt_edge*step # Determina o tamanho da aresta cúbica

    # ====================== Definições quanto ao output para o caso de um arquivo de dados txt ou resultado gráfico ==================== #
    if tp_arq == "txt":
        arq = open("{}/{}_{}.txt".format(path,tp,qnt), "w")
    elif tp_arq != "txt" and tp_arq != "":
        tipo = tp_arq
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
        X, Y, Z = raio*np.cos(u)*np.sin(v), raio*np.sin(u)*np.sin(v), raio*np.cos(v)
        # ======================================= Mantém a proporcionalidade na imagem final ============================================ #
        # ax.plot_surface(0.0001*np.cos(u)*np.sin(v) - step,0.0001*np.sin(u)*np.sin(v) - step,0.0001*np.cos(v) - step, color="red") # xyz=0
        ax.plot_surface(0.0001*np.cos(u)*np.sin(v) + edge,0.0001*np.sin(u)*np.sin(v),0.0001*np.cos(v), color="red") # x max
        ax.plot_surface(0.0001*np.cos(u)*np.sin(v),0.0001*np.sin(u)*np.sin(v) + edge,0.0001*np.cos(v), color="red") # y max
        ax.plot_surface(0.0001*np.cos(u)*np.sin(v),0.0001*np.sin(u)*np.sin(v),0.0001*np.cos(v) + edge, color="red") # z max
        # ============================================================= ## ============================================================== #
    # =============================================================== ## ================================================================ #

    qnt_part = 0 # Contagem das partículas no sistema

    # ============================================= Adiciona as partículas ao sistema =================================================== #
    while qnt_part < qnt:
        i,j,k = np.random.randint(qnt_edge), np.random.randint(qnt_edge), np.random.randint(qnt_edge)
        x,y,z = coordinates = [i*step, j*step, k*step]
        # coordinates = [x,y,z]
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
                        if not i%2 == 0:
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
        qnt_part += 1
        pos_w_part.append(coordinates)
        phi = np.random.rand()*2*np.pi
        theta = np.random.rand()*np.pi
        phi_u = np.random.rand()*2*np.pi
        theta_u = np.random.rand()*np.pi
        if tp_arq == "txt":
            if linha == 0: 
                pos = [f"{x} {y} {z} {phi} {theta}"]
                linha += 1
            else:
                pos = [f"\n{x} {y} {z} {phi} {theta}"]
            arq.writelines(pos)
        elif tp_arq != "txt" and tp_arq != "":
            ax.plot_surface(X + x, Y + y, Z + z,color="blue")
        if minimization:
            # for key in pos_values[:7]:
            #     pos_values[key].append(eval(key)) 
            pos_values["x"].append(x)
            pos_values["y"].append(y)
            pos_values["z"].append(z)
            pos_values["phi"].append(phi)
            pos_values["theta"].append(theta)
            pos_values["phi_u"].append(phi_u)
            pos_values["theta_u"].append(theta_u)
        perc = qnt_part/qnt*100
        print(f"Calculando as posições: [----- {perc:.0f}% -----]", end= "\r" if perc < 100 else "\n")
    # =============================================================== ## ================================================================ #

    # ========= Finalização do sistema inicial ========= #
    if tp_arq == "txt":
        arq.close()
    elif tp_arq != "txt" and tp_arq != "":
        print("Gerando imagem...")
        plt.savefig("{}/{}.{}".format(path,tp,tipo))
    # ======================= ## ======================= #

    # ====================================== Adiciona as contribuições energéticas em cada partícula ====================================== #
    if minimization:
        if tp_hys == "txt":
            graph = open(f"hyst_{tp}_{qnt}_{h_max}_{he_x}.txt", "w")
        else:
            fig, ax = plt.subplots(figsize=(10,10))
            # plt.axis("equal")
            ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.6)
            ax.set_ylim((-1, 1))
            ax.set_xlim((-h_max, h_max))
            ax.set_xlabel('h', size=14)
            ax.set_ylabel('m', size=14)
            ax.set_title(f"Curva de Histerese: {tp} {qnt} partículas")
        # Theta, Phi = sp.symbols("Theta Phi")
        # E = 0
        # m = 1
        # H = 3
        # He = 1
        # H = [He*k for k in [0,0,1]]
        # MM = [m*sp.sin(Theta)*sp.cos(Phi), m*sp.sin(Theta)*sp.sin(Phi), m*sp.cos(Theta)]
        # Ea = (((sp.sin(Theta)**2)*(sp.sin(2*Phi)**2))/4 + sp.cos(Theta)**2)*sp.sin(Theta)**2
        # Ez = -np.dot(MM,H)
        # Ed = 0
        p = 0.05
        L = int(h_max/p)
        # H = [0]*(5*L+1)
        H = [0]
        m = 1.0
        K = 1.0
        b = 0.1
        matrix_mag = []
        matrix_h = []
        start = time.time()
        # for i in range(1,L+1):
        #     H[i]=H[i-1] + p
        # for i in range(L+1,3*L+1):
        #     H[i]=H[i-1] - p 
        # for i in range(3*L+1,5*L+1):
        #     H[i]=H[i-1] + p
        for i in range(1,5*L+1):
            P = -p if 3*L+1 > i >= L+1 else p
            H.append(H[-1] + P)
        n = 0
        print(f"Calculando a contribuição energética das partículas: [----- 0% -----]", end= "\r")
        for h in H:
            for part in range(qnt):
                particle = [pos_values["x"][part], pos_values["y"][part], pos_values["z"][part]]
                theta = pos_values["theta"][part]
                phi = pos_values["phi"][part]
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
                            if particle != test_particle and test_particle in pos_w_part:
                                dist = np.sqrt((particle[0] - x)**2 + (particle[1] - y)**2 + (particle[2] - z)**2)
                                if dist < 2.1*raio:
                                    hi_x += he_x*np.sin(theta)*np.cos(phi)
                                    hi_y += he_x*np.sin(theta)*np.sin(phi)
                                    hi_z += he_x*np.cos(theta)
                # for key in pos_values[-3:]:
                #     pos_values[key][part] = eval(key)
                pos_values["hi_x"][part] = hi_x
                pos_values["hi_y"][part] = hi_y
                pos_values["hi_z"][part] = hi_z

            for part in range(qnt):
                # phi, theta, fu, tu, hi_x,hi_y,hi_z = [pos_values[key][part] for key in pos_values[3:]]
                tu = pos_values["theta_u"][part]
                fu = pos_values["phi_u"][part]
                theta = pos_values["theta"][part]
                phi = pos_values["phi"][part]
                hi_x = pos_values["hi_x"][part]
                hi_y = pos_values["hi_y"][part]
                hi_z = pos_values["hi_z"][part]

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
                    theta = theta - b*dgt
                    phi = phi - b*dgf
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
                pos_values["theta"][part] = theta
                pos_values["phi"][part] = phi
                n += 1
            s1 = 0.0
            s2 = 0.0
            for part in range(qnt):
                s1 += np.cos(pos_values["theta"][part])*np.sin(pos_values["theta_u"][part])
                s2 += np.sin(pos_values["theta_u"][part])
            mag = s1/s2
            # mag_txt.write(f"{i} {mag}\n")
            matrix_h.append(h)
            matrix_mag.append(mag)
            perc = n/(len(H)*qnt)*100
            print(f"Calculando a contribuição energética das partículas: [----- {perc:.0f}% -----]", end= "\r" if perc < 100 else "\n")
            # for k in range(11): # acima, abaixo e na mesma altura da partícula escolhida
            #     z = particle[2] + step*(k - 5)
            #     for j in range(11):
            #         y = particle[1] + step*(j - 5)
            #         for i in range(11):
            #             x = particle[0] + step*(i - 5)
            #             test_particle = [x,y,z]
            #             if particle[:3] != test_particle and test_particle in pos_w_part:
            #                 part_2 = pos_w_part.index(test_particle)
            #                 particle_2 = [x, y, z, pos_values["phi"][part_2], pos_values["theta"][part_2]]
            #                 dist = np.sqrt((particle[0] - x)**2 + (particle[1] - y)**2 + (particle[2] - z)**2)
            #                 if dist < 2.05*raio:
            #                     E_ex = 0
            #             else:
            #                 continue
            #perc = part/qnt*100
            #print(f"Calculando a contribuição energética das partículas: [----- {perc:.0f}% -----]", end= "\r" if perc < 100 else "\n")

        # ========================= Plot da curva de histerese ========================= #
        if tp_hys == "txt":
            for el in range(len(matrix_h)):
                if el == 0:
                    graph.write(f"{matrix_h[el]} {matrix_mag[el]}")
                else:
                    graph.write(f"\n{matrix_h[el]} {matrix_mag[el]}")
            graph.close()
        else:
            plt.plot(matrix_h,matrix_mag,color="red",linewidth="1")
            fig.savefig(f"{path}/{tp}/hys_{tp}_{qnt}_{h_max}_{he_x:.1f}.{tp_hys}",bbox_inches="tight")
        # ===================================== ## ===================================== #

    # ================================================================ ## ================================================================= #
    
    cube_vol = edge**3 
    part_vol = qnt*4*np.pi*raio**3/3
    print(f"Concentração: {(part_vol/cube_vol)*100:.2f}%")
    # mag_txt.close()
    end = time.time()
    print(f"Tempo: {round(end - start)} s")

for i in range(4,6):
    cubic_structure("BCC",3000,"",1,True,2.2,-0.2*i,"pdf")
import os, time
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

def magnetic_system(structure_type, num_of_particles, structure_file_format, concentration_of_particles,
                    minimization, max_value_mag_field, intensity_exchange, hys_file_format):

    radius = 1
    path = os.path.dirname(os.path.realpath(__file__))
    pos_w_part = []
    pos_values = {
                    "x": [],
                    "y": [],
                    "z": [], 
                    "phi": [], 
                    "theta": [], 
                    "phi_u": [], 
                    "theta_u": [], 
                    "hi_x": [], 
                    "hi_y": [], 
                    "hi_z": []
                }
    pos_values_keys = [key for key in pos_values]
    aprox = 0.09 
    start = time.time()
    # ====== Define a separação entre as partículas de acordo com o modelo cúbico ====== #
    if structure_type == "SC":
        step = 2*radius
        qnt_edge = ((2 - concentration_of_particles)*num_of_particles)**(1/3) 
    elif structure_type == "FCC":
        step = radius*np.sqrt(2)
        qnt_edge = (2*(2 - concentration_of_particles)*num_of_particles)**(1/3)
    elif structure_type == "BCC":
        step = 2*radius/np.sqrt(3)
        qnt_edge = (4*(2 - concentration_of_particles)*num_of_particles)**(1/3)
    qnt_edge = int(qnt_edge) if qnt_edge - int(qnt_edge) <= aprox else np.ceil(qnt_edge)
    edge = qnt_edge*step
    # ======================================= ## ======================================== #

    # ====================== Definições quanto ao output para o caso de um arquivo de dados txt ou resultado gráfico ==================== #
    if structure_file_format == "txt":
        arq = open(f"{path}/{structure_type}_{num_of_particles}.txt", "w")
    elif structure_file_format != "txt" and structure_file_format != "":
        tipo = structure_file_format
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
        X, Y, Z = radius*np.cos(u)*np.sin(v), radius*np.sin(u)*np.sin(v), radius*np.cos(v)

    # =============================================================== ## ================================================================ #

    # ============================================= Adiciona as partículas ao sistema =================================================== #
    qnt_part = 0 # Contagem das partículas no sistema
    while qnt_part < num_of_particles:
        i,j,k = np.random.randint(qnt_edge), np.random.randint(qnt_edge), np.random.randint(qnt_edge)
        x,y,z = coordinates = [i*step, j*step, k*step]

        # ====== Organiza as partículas de acordo com o modelo escolhido ====== #
        if coordinates in pos_w_part:
            continue
        else:
            if structure_type == "SC":
                pass
            elif structure_type == "FCC":
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
            elif structure_type == "BCC":
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
        # ================================= ## ================================ #

        qnt_part += 1

        # ======== Adiciona as coordenadas e orientação das partículas ======== #
        pos_w_part.append(coordinates)
        phi = np.random.rand()*2*np.pi
        theta = np.random.rand()*np.pi
        phi_u = np.random.rand()*2*np.pi
        theta_u = np.random.rand()*np.pi
        hi_x = 0
        hi_y = 0
        hi_z = 0
        # ================================= ## ================================ #

        # ============ Adiciona os elementos no arquivo de saída ============ #
        if structure_file_format == "txt":
            if len(pos_values["x"]) == 0: 
                pos = [f"{x} {y} {z} {phi} {theta}"]
            else:
                pos = [f"\n{x} {y} {z} {phi} {theta}"]
            arq.writelines(pos)
        elif structure_file_format != "txt" and structure_file_format != "":
            ax.plot_surface(X + x, Y + y, Z + z,color="blue")
        # ================================ ## =============================== #

        # ==== Se true adiciona as coordenadas às listas do dicionário pos_values ==== #
        if minimization:
            for key in pos_values:
                pos_values[key].append(eval(key))
        # ==================================== ## ==================================== #

        # ======================= Retorno do progresso do cálculo das posições ======================= #
        perc = qnt_part/num_of_particles*100
        print(f"Calculando as posições: [----- {perc:.0f}% -----]", end= "\r" if perc < 100 else "\n")
        # ============================================ ## ============================================ #

    # =============================================================== ## ================================================================ #

    # ============================== Finalização do sistema inicial ============================== #
    if structure_file_format == "txt":
        arq.close()
    elif structure_file_format != "txt" and structure_file_format != "":
        print("Gerando imagem...")
        plt.savefig(f"{path}/{structure_type}/system_{num_of_particles}.{structure_file_format}")
    # ============================================ ## ============================================ #

    # ====================================== Adiciona as contribuições energéticas em cada partícula ====================================== #
    if minimization:

        # =================================== Definições iniciais do output da curva de histerese =================================== #
        if hys_file_format == "txt":
            graph = open(f"hyst_{structure_type}_{num_of_particles}_{max_value_mag_field}_{intensity_exchange}.txt", "w")
        else:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.6)
            ax.set_ylim((-1, 1))
            ax.set_xlim((-max_value_mag_field, max_value_mag_field))
            ax.set_xlabel('h', size=14)
            ax.set_ylabel('m', size=14)
            ax.set_title(f"Curva de Histerese: {structure_type} {num_of_particles} partículas")
        # ============================================================ ## =========================================================== #

        # Theta, Phi, He = sp.symbols("Theta Phi He")
        # E = 0
        # m = 1
        # H = [He*k for k in [0,0,1]] # Vetor campo magnético externo
        # magnetic_moment = [m*sp.sin(Theta)*sp.cos(Phi), m*sp.sin(Theta)*sp.sin(Phi), m*sp.cos(Theta)]
        # Ea = (((sp.sin(Theta)**2)*(sp.sin(2*Phi)**2))/4 + sp.cos(Theta)**2)*sp.sin(Theta)**2
        # Zeemann_energy = -np.dot(MM,H)
        # Ed = 0
        step_magnetic_field = 0.05
        L = int(max_value_mag_field/step_magnetic_field)
        H = [0]
        m = 1.0
        K = 1.0
        b = 0.1
        matrix_mag = []
        matrix_h = []
        for i in range(1,5*L+1):
            S = -step_magnetic_field if 3*L+1 > i >= L+1 else step_magnetic_field
            H.append(H[-1] + S)
        
        # ======================= Cálculo e minimização da energia associada a cada partícula ======================= #
        n = 0 # Número de interações
        print(f"Calculando a contribuição energética das partículas: [----- 0% -----]", end= "\r")
        for h in H:
            for part in range(num_of_particles):
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
                                if dist < 2.1*radius:
                                    hi_x += intensity_exchange*np.sin(theta)*np.cos(phi)
                                    hi_y += intensity_exchange*np.sin(theta)*np.sin(phi)
                                    hi_z += intensity_exchange*np.cos(theta)
                for key in pos_values_keys[-3:]:
                    pos_values[key][part] = eval(key)

            for part in range(num_of_particles):
                phi, theta, fu, tu, hi_x,hi_y,hi_z = [pos_values[key][part] for key in pos_values_keys[3:]]

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
            s1 = 0
            s2 = 0
            for part in range(num_of_particles):
                s1 += np.cos(pos_values["theta"][part])*np.sin(pos_values["theta_u"][part])
                s2 += np.sin(pos_values["theta_u"][part])
            mag = s1/s2

            matrix_h.append(h)
            matrix_mag.append(mag)

            perc = n/(len(H)*num_of_particles)*100
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
        if hys_file_format == "txt":
            for el in range(len(matrix_h)):
                if el == 0:
                    graph.write(f"{matrix_h[el]} {matrix_mag[el]}")
                else:
                    graph.write(f"\n{matrix_h[el]} {matrix_mag[el]}")
            graph.close()
        else:
            plt.plot(matrix_h,matrix_mag,color="red",linewidth="1")
            fig.savefig(f"{path}/{structure_type}/hys_{num_of_particles}_{max_value_mag_field}_{intensity_exchange:.1f}.{hys_file_format}",bbox_inches="tight")
        # ===================================== ## ===================================== #

    # ================================================================ ## ================================================================= #
    
    cube_vol = edge**3 
    part_vol = num_of_particles*4*np.pi*radius**3/3
    print(f"Concentração: {(part_vol/cube_vol)*100:.2f}%")
    end = time.time()
    print(f"Tempo: {round(end - start)} s")

magnetic_system("SC",500,"png",1,True,3,0.4,"pdf")
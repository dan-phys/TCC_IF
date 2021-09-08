from structures import cubic_structure as cs
import sympy as sp
import numpy as np
import os, time

path = os.path.dirname(os.path.realpath(__file__))
txt = open("{}/teste.txt".format(path), "w")
values = {"x": [], "y": [], "z": [], "phi": [], "theta": []}
qnt = 1000

start = time.time()
cs("BCC",1,qnt,"",0.4,values['x'],values['y'],values['z'],values['phi'],values['theta'])

Theta, Phi = sp.symbols("Theta Phi")
f = (((sp.sin(Theta)**2)*(sp.sin(2*Phi)**2))/4 + sp.cos(Theta)**2)*sp.sin(Phi)**2
df_Theta = sp.diff(f,Theta) # Derivada parcial em relação a theta
df_Phi = sp.diff(f,Phi) # Derivada parcial em relação a phi
line = 0

for part_num1 in range(len(values['x'])):
    x_1 = values['x'][part_num1], y_1 = values['y'][part_num1], z_1 = values['z'][part_num1]
    for part_num2 in range(part_num1 + 1, len(values['x'])):
        x_2 = values['x'][part_num2], y_2 = values['y'][part_num2], z_2 = values['z'][part_num2]
        dist = np.sqrt()
    # text = [f"{round(values['x'][part_num],2)}", " ", f"{round(values['y'][part_num],2)}", " ", f"{round(values['z'][part_num],2)}"]
    # if line == 0:
    #     txt.writelines(text)
    #     line += 1
    # else:
    #     txt.writelines(["\n"] + text)
    # perc = part_num/qnt*100
    # print(f"Análise: [----- {perc:.2f}% -----]", end= "\r" if perc < 100 else "\n")

# print(f"Minimização: [----- {100.00}% -----]")
end = time.time()
txt.close()
print(end-start)
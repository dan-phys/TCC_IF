from structures import cubic_structure as cs
from sympy import *
import numpy as np
import os, time

path = os.path.dirname(os.path.realpath(__file__))
txt = open("{}/teste.txt".format(path), "w")
values = {"x": [], "y": [], "z": [], "phi": [], "theta": []}
qnt = 1000

start = time.time()
cs("BCC",1,qnt,"",0.4,values['x'],values['y'],values['z'],values['phi'],values['theta'])

Theta, Phi = symbols("Theta Phi")
E_a = (((sin(Theta)**2)*(sin(2*Phi)**2))/4 + cos(Theta)**2)*sin(Phi)**2
E = E_a
dE_Theta = diff(E,Theta) # Derivada parcial em relação a theta
dE_Phi = diff(E,Phi) # Derivada parcial em relação a phi
line = 0

for part_num in range(len(values['x'])):
    dET = dE_Theta.subs([(Theta, values['theta'][part_num]), (Phi, values['phi'][part_num])]).evalf()
    dEP = dE_Phi.subs([(Theta, values['theta'][part_num]), (Phi, values['phi'][part_num])]).evalf()
    # text = [f"{round(values['x'][part_num],2)}", " ", f"{round(values['y'][part_num],2)}", " ", f"{round(values['z'][part_num],2)}"]
    # if line == 0:
    #     txt.writelines(text)
    #     line += 1
    # else:
    #     txt.writelines(["\n"] + text)
    # perc = part_num/qnt*100
    # print(f"Análise: [----- {perc:.2f}% -----]", end= "\r" if perc < 100 else "\n")

end = time.time()
txt.close()
print(end-start)
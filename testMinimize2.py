import numpy as np

angles = (80*np.pi/180,15*np.pi/180)
a = 0.5
n = 0
x,y = angles
mx = np.sin(x)*np.cos(y)
my = np.sin(x)*np.sin(y)
mz = np.cos(x)        
dmxt = np.cos(x)*np.cos(y)
dmxf = -np.sin(x)*np.sin(y)
dmyt = np.cos(x)*np.sin(y)
dmyf = np.sin(x)*np.cos(y)
dmzt = -np.sin(x)
dmzf = 0.00000000000
deut=-(mx**3*dmxt + my**3*dmyt + mz**3*dmzt)
deuf=-(mx**3*dmxf + my**3*dmyf + mz**3*dmzf)
gt =  deut
gf = deuf
g = np.sqrt((gt)**2+(gf)**2)
while g >= 10**(-4):
    x = x - gt*a
    y = y - gf*a
    mx = np.sin(x)*np.cos(y)
    my = np.sin(x)*np.sin(y)
    mz = np.cos(x)        
    dmxt = np.cos(x)*np.cos(y)
    dmxf = -np.sin(x)*np.sin(y)
    dmyt = np.cos(x)*np.sin(y)
    dmyf = np.sin(x)*np.cos(y)
    dmzt = -np.sin(x)
    dmzf = 0.00000000000
    deut=-(mx**3*dmxt + my**3*dmyt + mz**3*dmzt)
    deuf=-(mx**3*dmxf + my**3*dmyf + mz**3*dmzf)
    gt =  deut
    gf = deuf
    g = np.sqrt((gt)**2+(gf)**2)
    n += 1
print(x, y)
print(n)
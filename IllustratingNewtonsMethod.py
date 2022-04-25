from numpy import *
%matplotlib inline
import matplotlib.pyplot as plt
nr = 1000
nc = 1000
r = 1.5
x = linspace(-r,r,nc)
y = linspace(r,-r,nr)
x;
X,Y = meshgrid(x,y)
Z = X + 1j*Y
def g(z): return 2/3*z + 1/3/z**2
a = zeros((nr,nc,3),dtype = uint8)

red =    a[:,:,0]
blue =   a[:,:,1] 
green =  a[:,:,2]

redroot   =  1.
blueroot  = -1/2 + 1j*sqrt(3)/2
greenroot = -1/2 - 1j*sqrt(3)/2

for k in range(10):
    
    red   [ abs(Z-redroot)   < 0.1 ] += 1  #add one to every pixel where Z now close to red root
    blue  [ abs(Z-blueroot)  < 0.1 ] += 1  #add one to every pixel where Z now close to blue root
    green [ abs(Z-greenroot) < 0.1 ] += 1  #add one to every pixel where Z now close to green root
    
    Z = g(Z) # i just wanna replace Z by g(Z)

from scipy.misc import imsave
a = array(a/a.max()*255,dtype = uint8)
plt.figure(figsize = (10,10))
plt.imshow(a,interpolation = 'none',extent = [-r,r,-r,r]);
plt.figure(facecolor = 'w')
plt.subplot(111,aspect = 1)

t = linspace(0,2*pi,60)
c = cos(t)
d = sin(t)
#plt.plot(c,d,'k',alpha = 0.25)

r = 0.8
c *= r
d *= r
z = c + d*1j

w = z**3-1

plt.plot(z.real,z.imag,'ro', alpha = 0.5)
plt.plot(w.real, w.imag,'co')
for e,f in zip(z,w):
    plt.plot([e.real,f.real],[e.imag,f.imag],'g',alpha=0.25)
plt.figure(facecolor = 'w')
plt.subplot(111,aspect = 1)

t = linspace(0,2*pi,60)
c = cos(t)
d = sin(t)
#plt.plot(c,d,'k',alpha = 0.25)

r = 0.8
c *= r
d *= r
z = c + d*1j

w = 2/3*z + 1/3/z**2

plt.plot(z.real,z.imag,'ro', alpha = 0.5)
plt.plot(w.real, w.imag,'co')
for e,f in zip(z,w):
    plt.plot([e.real,f.real],[e.imag,f.imag],'g',alpha=0.25)
j = []
k = []
def newtonmethod(f,fprime,z,tol): 
    while True:
        print(z)
        j.append(z)
        s = - f(z)/fprime(z)
        k.append(z+s)
        z += s
        if abs(s)<tol: return z
        for h,i in zip(j,k):
            plt.plot([h.real,i.real],[h.imag,i.imag],'r',ms = 5, alpha = 1)
def f(z): return z**3 - 1
def fprime(z): return 3*(z**2)
z=0
newtonmethod(f,fprime,0.5 + 0.5*1j,1.e-12)

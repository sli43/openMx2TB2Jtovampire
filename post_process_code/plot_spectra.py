import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors

plt.rcParams["font.family"]="Arial"
plt.rcParams["font.size"]=16
Nx=16
Ny=16
Nz=16
data=np.transpose(np.loadtxt('spectral.dat',unpack=True))

size1=data.shape[0]
size2=data.shape[1]-1
print(size1, size2)

x=np.zeros([size1,size2])
y=np.zeros([size1,size2])
z=np.zeros([size1,size2])
maxpoint=np.zeros([size2,3])

for j in range(1,size2+1):
   f=data[:,j]
   maxpoint[j-1,0]=data[np.argmax(f),0]
   maxpoint[j-1,1]=j-1
   maxpoint[j-1,2]=np.max(f)

for i in range(0,size1):
 for j in range(1,size2+1):
     x[i,j-1]=data[i,0]
     y[i,j-1]=(j-1)

     z[i,j-1]=data[i,j]/maxpoint[j-1,2]

fig=plt.figure(figsize=(6,6))
ax1=plt.subplot(111)
ax1.set_position([0.15,0.15,0.7,0.6])

plt.pcolormesh(np.transpose(y),np.transpose(x),np.transpose(z),cmap='hot',shading='gouraud',vmin=0,vmax=1)
plt.plot(maxpoint[:,1],maxpoint[:,0],'sc',ms=3)
plt.plot([Nx/2,Nx/2],[0,1000],'--w',lw=2)
plt.plot([(Nx+Ny)/2,(Nx+Ny)/2],[0,1000],'--w',lw=2)
plt.plot([(Nx+Ny+Nz)/2,(Nx+Ny+Nz)/2],[0,1000],'--w',lw=2)

plt.ylim(0,800)
plt.xticks([0,Nx/2,(Nx+Ny)/2,(Nx+Ny+Nz)/2,(Nx+Ny+2*Nz)/2],['G','X','M','R','G'])

plt.ylabel('$E(q)$ [meV]')
fig.savefig('T=0.1.pdf')
plt.show()

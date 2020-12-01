#environment, only 'matplotlib' and 'numpy' are needed
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from matplotlib import animation
import random
import os


N=0 #total number of particles
hashTable=[] #hashTable, array(not actived now)
x=[] #cordination x, array
y=[] #cordination y
vx=[] #velocity x
vy=[] #velocity y
p=[] #pressure, array
rho=[] #density
m=[] #mass of particles
bx=0 #boundary in x-axis, 0-bx
by=0 #0-by
k=0 #size of a hash square (not actived)
hx=0 #number of hashTable's square in x-axis 
hy=0 #number of hashTable's square in y-axis 
tau=0 #time step
s=0 #number of boundary nodes

### not used variables
rho0=0 #reference density
gamma=0 #used in p-rho equation
B=0 #used in p-rho equation


def appx(xi,yi,temp,r,k,t): #generate horizonal boundary
    xi=xi/t
    yi=yi/t
    for j in range(0,temp,2):
        x.append(xi+0.01+j/t)
        y.append(yi+0.5)
        m.append(1)
        vx.append(0)
        vy.append(0)
        rho.append(r)
        p.append(k)

        x.append(xi+0.01+j/t)
        y.append(yi+1)
        m.append(1)
        vx.append(0)
        vy.append(0)
        rho.append(r)
        p.append(k)

        x.append(xi+0.01+j/t)
        y.append(yi+0.02)
        m.append(1)
        vx.append(0)
        vy.append(0)
        rho.append(r)
        p.append(k)

    return x,y,m,vx,vy,rho,p

def appy(xi,yi,temp,r,k,t):  #generate verticle boundary
    xi=xi/t
    yi=yi/t
    for j in range(0,temp,2):
        x.append(xi+0.5)
        y.append(yi+0.01+j/t)
        m.append(1)
        vx.append(0)
        vy.append(0)
        rho.append(r)
        p.append(k)

        x.append(xi+0.015)
        y.append(yi+0.01+j/t)
        m.append(1)
        vx.append(0)
        vy.append(0)
        rho.append(r)
        p.append(k)

        x.append(xi+1)
        y.append(yi+0.01+j/t)
        m.append(1)
        vx.append(0)
        vy.append(0)
        rho.append(r)
        p.append(k)
    return x,y,m,vx,vy,rho,p

def inilization0():  #generate all particles
    bx=100
    by=100
    k=0 #initial pressure
    r=1 #initial density
    temp=int(bx/3)-1 #temporary various
    t=10 #control the size of domain(placed on bx&by, no need to pay attention on)

    ### generate boundary particles
    x,y,m,vx,vy,rho,p=appx(5,by,bx,r,k,t) 
    x,y,m,vx,vy,rho,p=appx(5,0,bx,r,k,t)
    x,y,m,vx,vy,rho,p=appy(0,5,by,r,k,t)
    x,y,m,vx,vy,rho,p=appy(bx,5,by,r,k,t)
    s=len(x) #record the number of boundary particles

    ### generate moving particles
    for i in range(temp,temp*2,4):
        for j in range(temp,temp*2,4): 
            ### add one partile and all its information ###
            x.append(i/t+0.5+random.random()/5) 
            y.append(j/t+0.5/t+0.5+random.random()/5) 
            m.append(1) 
            vx.append(0) 
            vy.append(0) 
            rho.append(1) 
            p.append(k) 
    # x.append(2)
    # y.append(3)
    # m.append(1)
    # vx.append(0)
    # vy.append(0)
    # rho.append(1)
    # p.append(0)

    # x.append(9)
    # y.append(9)
    # m.append(1)
    # vx.append(0)
    # vy.append(0)
    # rho.append(1)
    # p.append(0)

    N=len(x)
    return x,y,vx,vy,rho,p,N,bx,by,s #return the valuables

def inilization(hx,hy): #inilize the 2-D hashTahle (not actived)
    temp=()
    hashTable=[]
    for i in range(hx):
        hashTable.append([])
        for j in range(hy):
            hashTable[i].append(temp)
    return hashTable

def hashSize(bx,by,k): #set the sizes of the hashTahle (not actived)
    hx=int(bx/k)+4  # domain+2 additional space
    hy=int(by/k)+4
    n=hx*hy
    return hx,hy,n

def hashIn(i,hashTable): # insert particles into corresponding hash square (not actived)
    ix,iy=returnHashIndex(i)
    temp=(i,) #( Each element in hashTable is a list recording the index of particles)
    hashTable[ix][iy]=hashTable[ix][iy]+temp
    return hashTable

def returnHashIndex(i): #return the address of a particle in hashTable
    t=10
    ix=x[i]*t
    iy=y[i]*t
    ix=int(ix/k)+1
    iy=int(iy/k)+1
    return ix, iy

def kernalFunction(i,j): #kernalFunction
    temp=((x[i]-x[j])**2+(y[i]-y[j])**2)**0.5
    w=0
    if ((x[i]-x[j])**2+(y[i]-y[j])**2>4):
        w=0
    elif ((x[i]-x[j])**2+(y[i]-y[j])**2>1):
        w=5/14/np.pi*(2-temp)**3
    elif (temp>=0):
        w=5/14/np.pi*(2-temp)**3-5/14/np.pi*(1-temp)**3*4
    return w

def kernalFunctionX(i,j): #particle difference of kernalFunction in x
    temp=((x[i]-x[j])**2+(y[i]-y[j])**2)**0.5 #i, j: the index of particle i, j
    w=0
    if (temp>2):
        w=0
    elif (temp>1):
        w=15/14/np.pi/temp*(x[i]-x[j])*(2-temp)**2
    elif (temp>0):
        w=15/14/np.pi/temp*(x[i]-x[j])*(2-temp)**2-15/14/np.pi/temp*(x[i]-x[j])*(1-temp)**2*4
    return w

def kernalFunctionY(i,j): #particle difference of kernalFunction in y
    temp=((x[i]-x[j])**2+(y[i]-y[j])**2)**0.5
    w=0
    if (temp>2):
        w=0
    elif (temp>1):
        w=15/14/np.pi/temp*(y[i]-y[j])*(2-temp)**2
    elif (temp>0):
        w=15/14/np.pi/temp*(y[i]-y[j])*(2-temp)**2-15/14/np.pi/temp*(y[i]-y[j])*(1-temp)**2*4
    return w

def kernalSearch(hashTable, i): #return a list of surrounding particles of particle i
    x, y=returnHashIndex(i)
    neighbour=hashTable[x][y]+hashTable[x-1][y]+hashTable[x][y-1]
    neighbour=neighbour+hashTable[x-hx+1][y]+hashTable[x][y-hy+1]
    neighbour=neighbour+hashTable[x-hx+1][y-1]+hashTable[x-1][y-hy+1]
    neighbour=neighbour+hashTable[x-hx+1][y-hy+1]+hashTable[x-1][y-1]
    return neighbour

def updateDensityPressure(rho,p,s): #update Density & Pressure of all particles
    rho_hat=rho.copy() #make a copy of density
    p_hat=p.copy()
    for i in range(N): 
        # neighbour=kernalSearch(hashTable, i)
        rho_temp=0   #inilization of temporary variable
        # p_hat[i]=B*((rho[i]/rho0)**gamma-1)
        p_hat[i]=rho[i]**4 # density-pressure function
        # for j in neighbour:
        for j in range(N):
            # rho_temp=rho_temp+rho[i]*m[j]/rho[j]*(vx[i]-vx[j])*kernalFunctionX(i,j)+rho[i]*m[j]/rho[j]*(vy[i]-vy[j])*kernalFunctionY(i,j)
            rho_temp=rho_temp+m[j]*kernalFunction(i,j)  # demsity equation (accumalate)
        rho_hat[i]=rho_temp
    return rho_hat, p_hat

def updatePosition(x,y,vx,vy,tau): #update cordinate
    for i in range(s,N):
        x[i]=x[i]+vx[i]*tau
        y[i]=y[i]+vy[i]*tau
    return x, y

def updateVelocity(x,y,vx,vy,s): #update velocity
    vx_hat=vx.copy()
    vy_hat=vy.copy()
    for i in range(s,N): #skip the boundary particles
        # neighbour=kernalSearch(hashTable, i)
        vx_temp=0
        vy_temp=0
        # for j in neighbour:
        for j in range(N):
            vx_temp=vx_temp+m[j]*(p[i]/rho[i]/rho[i]+p[j]/rho[j]/rho[j])*kernalFunctionX(i,j)
            vy_temp=vy_temp+m[j]*(p[i]/rho[i]/rho[i]+p[j]/rho[j]/rho[j])*kernalFunctionY(i,j)
        vx_hat[i]=vx_temp*tau+vx_hat[i]
        vy_hat[i]=vy_temp*tau+vy_hat[i]
    return vx_hat, vy_hat




### inilization of parament
x,y,vx,vy,rho,p,N,bx,by,s=inilization0()
gamma=10
rho0=5
tau=0.01
B=7
k=5
hx,hy,n=hashSize(bx,by,k)

### inilization of HashTable (not actived)
hashTable=inilization(hx,hy)
for i in range(N):
    hashTable=hashIn(i,hashTable)

###plot the figure
plt.figure()
plt.scatter(x[s:],y[s:])
plt.scatter(x[:s],y[:s])
print(s)

### itelation times
ite=20

## pre-work
rho,p=updateDensityPressure(rho,p,s)
vx_hat,vy_hat=updateVelocity(x,y,vx,vy,s)
for i in range(N):
    vx[i]=(vx_hat[i]+vx[i])/2
    vy[i]=(vy_hat[i]+vy[i])/2


### itelation (ATTENTION, figures are saved in ./pn folder)
for i in range(ite):
    x,y=updatePosition(x,y,vx,vy,tau)
    rho,p=updateDensityPressure(rho,p,s)
    vx,vy=updateVelocity(x,y,vx,vy,s)
    # print(i,vx[s],vy[s],rho[s],p[s])

    ### plot & save all the figures 
    path='./pn/fig'+str(i)+'.jpg'
    plt.figure()
    plt.scatter(x[s:],y[s:])
    plt.scatter(x[:s],y[:s])
    plt.savefig(path)


    ### Used for update hashTable (not actived)
    # hashTable=inilization(hx,hy)
    # for j in range(N):
    #     hashTable=hashIn(j,hashTable)

















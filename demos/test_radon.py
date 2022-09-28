#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import io
import numpy as np
import pyradon as pr
#import scipy
# import radon_op
import matplotlib.pyplot as plt
# import pcg
# 
# REFERENCE
# Chen, 2018, GEO, Automatic velocity analysis using high-resolution hyperbolic Radon transform
# Zhang, 2022, GJI, Improving receiver function with high-resolution Radon transform


d = np.loadtxt('./data/test_z_data.txt')
h = np.loadtxt('./data/test_dist_data.txt')

nh = len(h)
dq = 0.005
q = np.arange(-2,2,dq)
nq = len(q)

t = np.linspace(-5.0,144.9,1500)
dt = 0.1
nt = len(t)
#m = np.zeros((nt,nq))

Param = {}

class Param:
    h=h
    v=1./q
    q=q
    nt=nt
    dt=dt
    typ=1
ma=np.zeros((nt,nq))

nx = len(d)

N1 = 10  # CG Iterations (Internal loop);
N2 = 1   # Update of weights for the sparse solution: N1 = 1 LS;  N2 > 3 for High Res (Sparse) solution


## Adjoint
m = pr.radon(d,Param,-1)
print("Adjoint is done\n")
## Inversion
mi,misfit = pr.pcg(pr.radon,Param,d,np.zeros(ma.shape),niter_in=N1,niter_out=N2,verb=1)
print("Inversion is done\n")

di = pr.radon(mi,Param,1)

plt.figure(figsize=(14,16)) # Plot original data with Radon spectrum (d,m)
plt.subplot(1,2,1)
d_max=np.max(d[:,10])/10
for n in range(np.size(d,1)):
    x=d[:,n]
    plt.plot(x/d_max+h[n-1],t,'k')
ax=plt.gca()
plt.xlim(np.min(h),np.max(h))
plt.ylim(np.min(t),np.max(t))
ax.invert_yaxis()
plt.xlabel('Distance(km)',fontdict={'weight':'normal','size': 20})
plt.ylabel('Time (s)',fontdict={'weight':'normal','size': 20})
plt.tick_params(labelsize=20)
plt.title('RFs',fontdict={'weight':'normal','size': 20})

plt.subplot(1,2,2)
plt.imshow(m/np.max(m),cmap='RdBu',clim=(-0.2,0.2),extent=[q.min(),q.max(),t.max(),t.min()],interpolation='nearest',aspect='auto')
#plt.xlim(0,3200)
#plt.ylim(1500,0)
ax=plt.gca()
ax.yaxis.set_ticks_position('right')
plt.xlabel('Slowness (s/km)',fontdict={'weight':'normal','size': 20})
#plt.ylabel('Time (s)',fontdict={'weight':'normal','size': 20})
plt.title('Radon Spectrum',fontdict={'weight':'normal','size': 20})
plt.tick_params(labelsize=20)
#plt.colorbar(extend='both',orientation='horizontal')
plt.savefig('test.png',bbox_inches='tight',dpi=300)



plt.figure(figsize=(14,16)) # Plot Radon denoised data with Radon spectrum (di,mi)
plt.subplot(1,2,1)
di_max=np.max(di[:,10])/10
for n in range(np.size(di,1)):
    x=di[:,n]
    plt.plot(x/di_max+h[n-1],t,'k')
ax=plt.gca()
plt.xlim(np.min(h),np.max(h))
plt.ylim(np.min(t),np.max(t))
ax.invert_yaxis()
plt.xlabel('Distance(km)',fontdict={'weight':'normal','size': 20})
plt.ylabel('Time (s)',fontdict={'weight':'normal','size': 20})
plt.tick_params(labelsize=20)
plt.title('RFs',fontdict={'weight':'normal','size': 20})

plt.subplot(1,2,2)
plt.imshow(mi/np.max(mi),cmap='RdBu',clim=(-0.2,0.2),extent=[q.min(),q.max(),t.max(),t.min()],interpolation='nearest',aspect='auto')
#plt.xlim(0,3200)
#plt.ylim(1500,0)
ax=plt.gca()
ax.yaxis.set_ticks_position('right')
plt.xlabel('Slowness (s/km)',fontdict={'weight':'normal','size': 20})
#plt.ylabel('Time (s)',fontdict={'weight':'normal','size': 20})
plt.title('Radon Spectrum',fontdict={'weight':'normal','size': 20})
plt.tick_params(labelsize=20)
#plt.colorbar(extend='both',orientation='horizontal')
plt.savefig('test_i.png',bbox_inches='tight',dpi=300)

import numpy as np
import matplotlib.pyplot as plt
import pyradon as pr
import time

def ricker(f,dt,tlength=None):
    # ricker: Ricker wavelet of central frequency f.
    #
    # INPUT:
    # f : central freq. in Hz (f <<1/(2dt) )
    # dt: sampling interval in sec
    # tlength : the duration of wavelet in sec
    #
    # OUTPUT: 
    # w:  the Ricker wavelet
    # tw: time axis
    #
    # Example
    #
    #   [w,tw] = ricker(10,0.004,0.2);
    #   pyplot.plot(tw,w);
    import numpy as np
    
    if tlength!=None:
        nw=np.floor(tlength/dt)+1;
    else:
        nw=2.2/f/dt;
        nw=2*np.floor(nw/2)+1;
    nc=np.floor(nw/2);
    nw=int(nw)
    w =np.zeros(nw);
    
    k=np.arange(1,nw+1,1);
    alpha = (nc-k+1)*f*dt*np.pi;
    beta=np.power(alpha,2);
    w = (1.-beta*2)*np.exp(-beta);
    tw = -(nc+1-k)*dt;
    return w,tw

dp=0.05;dh=0.1;

h=np.linspace(0,6,int(6/dh+1));

nh=len(h);
p=np.linspace(-0.6,0.60,int(1.21/dp)+1);
nnp=len(p);

dt=0.2;nt=501;
m=np.zeros([nt,nnp]);
d=np.zeros([nt,nh]);
w,tw=ricker(0.5,dt);


tau=[100,200,250];
p0=[0,0.1,0.1];


m[tau[0]-1,12]=1;
m[tau[1]-1,18]=0.4;
m[tau[2]-1,18]=0.5;
# t = np.convolve(m[:,1],w,'same');
import scipy.signal
w=np.expand_dims(w,1);
m=scipy.signal.convolve2d(m,w,'same');

par={'h':h,'v':1/(p+0.0000000000000001),'nt':nt,'dt':dt,'typ':1}

## Python version
d=pr.radon(m,par,1);

# ## C version
tic = time.perf_counter()
### C-version
#Forward operator
d_c=pr.radonc(m,par,1);
#Adjoint operator
ma_c = pr.radonc(d_c,par,-1);
toc = time.perf_counter()
print(f"C version takes {toc - tic:0.4f} seconds");

## LS
tic = time.perf_counter()
mi_c=pr.radon_l2c(par,d,np.zeros(m.shape),niter=10,verb=1)[0]
toc = time.perf_counter()
print(f"C version LS takes {toc - tic:0.4f} seconds");

## High-resolution 
tic = time.perf_counter()
ml1_c=pr.radon_l1c(par,d,np.zeros(m.shape),niter_in=10,niter_out=3,verb=1)[0]
toc = time.perf_counter()
print(f"C version L1 takes {toc - tic:0.4f} seconds");

plt.figure(figsize=(8, 8));
plt.subplot(1,5,1)
plt.imshow(m,aspect='auto');plt.title('Model');
plt.subplot(1,5,2)
plt.imshow(d,aspect='auto');plt.title('Data');plt.gca().set_yticks([]);
plt.subplot(1,5,3)
plt.imshow(ma_c,aspect='auto');plt.title('Adjoint');plt.gca().set_yticks([]);
plt.subplot(1,5,4)
plt.imshow(mi_c,aspect='auto');plt.title('LS inversion');plt.gca().set_yticks([]);
plt.subplot(1,5,5)
plt.imshow(ml1_c,aspect='auto');plt.title('L1 inversion');plt.gca().set_yticks([]);
plt.show()









### inversion-based Radon transform including Least-squares Radon and High-resolution Radon

import numpy as np
from .pcg import pcg
from .radon_op import radon
from radoncfun import *

def radon_l2(par,d,m0,niter=10,verb=0):
	'''
	Least-squares Radon
	
	Written on Oct, 5, 2022
	'''
	mi,misfit=pcg(radon,par,d,m0,niter_in=niter,niter_out=1,verb=verb);
	
	return mi,misfit


def radon_l1(par,d,m0,niter_in=10,niter_out=3,verb=0):
	'''
	High-resolution Radon
	
	Written on Oct, 5, 2022
	'''
	mi,misfit=pcg(radon,par,d,m0,niter_in=niter_in,niter_out=3,verb=verb);

	return mi,misfit
	
	
def radon_l2c(par,d,m0,niter=10,verb=0):
	'''
	Least-squares Radon implemented in C
	
	Written on Oct, 5, 2022
	'''
	
	h = par['h'];
	v = par['v'];
	nt = par['nt'];
	dt = par['dt'];
	typ = par['typ'];

	nh = len(h);
	nv = len(v);
	
	d=np.float32(d.flatten(order='F'));
	m0=np.float32(m0.flatten(order='F'));
	v=np.float32(v);
	h=np.float32(h);
	
	tmp=radonc_inv(d, m0, v, h, typ, niter, 1, nt, nv, nh, dt, verb);
	mi=tmp[:nt*nv].reshape(nt,nv,order='F');
	misfit=tmp[nt*nv:];
	
	return mi,misfit


def radon_l1c(par,d,m0,niter_in=10,niter_out=3,verb=0):
	'''
	Least-squares Radon implemented in C
	
	Written on Oct, 5, 2022
	'''
	
	h = par['h'];
	v = par['v'];
	nt = par['nt'];
	dt = par['dt'];
	typ = par['typ'];

	nh = len(h);
	nv = len(v);

	d=np.float32(d.flatten(order='F'));
	m0=np.float32(m0.flatten(order='F'));
	v=np.float32(v);
	h=np.float32(h);
	
	tmp=radonc_inv(d, m0, v, h, typ, niter_in, niter_out, nt, nv, nh, dt, verb);
	mi=tmp[:nt*nv].reshape(nt,nv,order='F');
	misfit=tmp[nt*nv:];
	
	return mi,misfit







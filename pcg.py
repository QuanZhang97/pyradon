#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import math
import numpy as np

def pcg(oper,par,d,m0,niter_in=10,niter_out=1,verb=0):
# pcg: Precondioned CG for solving sparsity-promoting inverse problems
# min || d - Fm||_2^2 + \mu || m ||_1
#
# INPUT
# operator: forward operator F
# Param:    parameter struct of the forward operator
# d:        RHS of the inverse problem
# m0:       initial model estimation
# Niter_in: inner iteration NO
# Niter_out:outer iteration NO
# verb:     verbosity
# 
# OUTPUT
# m:        estimated model
# misfit:   misfit history
# 
# REFERENCE
# Chen, 2018, GEO, Automatic velocity analysis using high-resolution hyperbolic Radon transform

    u=m0;
    P=np.ones(u.shape); #P is an diagonal weighting operator, so I used capital P
    kc=0;
    mis=[];
    m=u;
    for l in range(niter_out):
        di=oper(P*u,par,1);
        r=d-di;
        
        g=oper(r,par,-1);
        g=g*P;
        s=g;
        gammam=np.sum(g*np.conj(g));
        for k in range(niter_in):
            q=oper(P*s,par,1);
            den=np.sum(q*np.conj(q));
            alpha=gammam/(den+1.e-8);
            u=u+alpha*s;
            r=r-alpha*q;
            mis.append(np.sum(r*np.conj(r)));
            g=oper(r,par,-1);
            g=g*P;
            gamma=np.sum(g*np.conj(g));
            beta=gamma/(gammam+1.e-7);
            gammam=gamma;
            s=g+beta*s;
            if verb:
                print("Outer iteration=%d, Inner iteration = %d, Current misfit=%0.5g \n"%(l+1,k+1,mis[kc]));
            kc=kc+1;
        m=P*u;
        P=np.abs(m/np.max(m.flatten()))+0.001;
    mis=np.array(mis);
    return m,mis


def cgdot(inpt):
# Dot product

    temp = np.dot(inpt,np.conj(inpt))
    out = np.sum(temp)
    return out

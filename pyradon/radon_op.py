#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def radon(inpt,Param,operator):
#Forward and Adjoint operators for Radon RT in the time domain.
#  IN   	in:   		intput data 
#       	Param:  		parameters combination
#		Param.h:		offset
#		Param.v:		velocity
#		Param.nt:		number of samples
#		Param.dt:		time interval
#		Param.type:    1: linear 2: parablic 3: hyperbolic
#
#	   	operator: 
## 			operator =  1 means impute is m(tau,v) and output is d(t,x) FORWARD  OP
##			operator = -1 means input is d(t,x) and output is m(tau,v)  ADJOINT  OP 
#      
#  OUT   out:  		output data
# 
#  Copyright (C) 2014 The University of Texas at Austin
#  Copyright (C) 2014 Yangkang Chen
#
#  Example:
#  test/test_radon_recon_hyper.m
#  test/test_radon_recon_linear.m
#  test/test_radon_demul.m
#  
#  Dot test example:
#  nt=500;
#  nv=100;
#  nh=100;
#  h=1:nh;
#  dt=1;
#  type=1;
#  
#  Param.h=h;
#  Param.nt=nt;
#  Param.dt=dt;
#  Param.v=linspace(-5,10,nv);
#  Param.type=type;
# 
#  forward
#  m1 = randn(nt,nv); 
#  [d1 ] = radon_op(m1,Param,1);
#  adjoint
#  d2 = randn(nt,nh); 
#  [m2 ] = radon_op(d2,Param,-1);
# 
# REFERENCE
# Chen, 2018, GEO, Automatic velocity analysis using high-resolution hyperbolic Radon transform



    h = Param.h
    v = Param.v
    nt = Param.nt
    dt = Param.dt
    typ = Param.typ
    
    nh = len(h)
    nv = len(v)
    
    if operator == -1:
        m = np.zeros((nt,nv),dtype='float64')
    if operator ==  1:
        d = np.zeros((nt,nh),dtype='float64')
    
    if operator == -1:
        d = inpt
    if operator ==  1:
        m = inpt
    
    hmax=np.max(abs(h))
    
    for itau in range(nt):
        #print(itau)
        for ih in range(nh):
            for iv in range(nv):
## This can also be replaced by Parabolic or linear integration 
                if typ == 1:
                    t = itau*dt + h[ih]/v[iv]	
                    it = int(np.floor(t/dt)+1)     
                elif typ == 2:
                    t = itau*dt + h[ih]*h[ih]*v[iv]/hmax/hmax  #curvature
                    it = int(np.floor(t/dt)+1)
            #if(it<=0) it=1;end
                elif typ == 3:
                    t = np.sqrt ((itau*dt)^2 + (h[ih]/v[iv])^2 )
                    it = int(np.floor(t/dt)+1)
                else:
                    t = np.sqrt ((itau*dt)^2 + (h[ih]/v[iv])^2 )
                    it = int(np.floor(t/dt)+1)
                #print(it),print(ih),print(itau),print(iv,'\n')
                if (it+1<=nt) & (it+1>0):  
                    if (operator==-1):
                        m[itau,iv]=m[itau,iv]+d[it,ih]
                    if (operator==1):
                        d[it,ih]=d[it,ih]+m[itau,iv]
    if operator == 1:
        out = d
    if operator ==-1:
        out = m

    return out


close all; clear; clc

addpath('../functions/')

L = 38.6;
N = 64;
symm = true;
dt = 0.1;

[x,~] = domain(L,N);
u0 = sin(2*pi*x/L);

u_guess = field2vector(u0,N,symm);
T_guess = 50;

[u_best,T_best] = search4PO(u_guess,T_guess,dt,L,N,symm); 

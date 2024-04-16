close all; clear; clc

addpath('../functions/')

L = 38.6;
N = 64;
symm = true;
T = 750;
dt = 0.1;

[x,~] = domain(L,N);
u0 = sin(2*pi*x/L);

v0 = field2vector(u0,N,symm);
[vT,~] = KSE_integrate(v0,T,dt,0,L,N,symm);

uT = vector2field(vT,N,symm);

figure
    plot(x,uT,'LineWidth',2)
    hold on; grid on
    plot(x,u0,'LineWidth',2)
    xlabel('x'); ylabel('u')
    legend('u(T)','u(0)')

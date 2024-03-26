%%% DESCRIPTION -----------------------------------------------------------
%   Jacobian of the KSE flow map: J=du(t)/du(0)


%%% INPUTS ----------------------------------------------------------------
%   u0      reference point of the Jacobian (column state vector)
%   t       integration time interval
%   eps     perturbation magnitude for finite difference derivatrives
%   dt      step size in time integrations
%   L       domain length
%   N       spatial resolution
%   symm	center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   J       Jacobian matrix


function J = Jacobian(u0,t,eps,dt,L,N,symm)
    n = length(u0);
    J = zeros(n,n);
    
    %%% to be completed...
end
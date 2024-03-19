%%% DESCRIPTION -----------------------------------------------------------
%   dealiasing a 1D periodic field in spectral state using the 2/3 rule


%%% INPUTS ----------------------------------------------------------------
%   u   vector of Fourier coefficients (vector of complex numbers)


%%% OUTPUTS ---------------------------------------------------------------
%   v   the same as 'u' with 1/3 of the highest frequency modes set to 0


function v = dealiase(u)
    N = length(u);
    Nd = round(N/3);
    
    v = u;
    v(Nd+2:N-Nd) = 0;
end
%%% DESCRIPTION -----------------------------------------------------------
%   1D periodic field corresponding to a minimal vector of state variables


%%% INPUTS ----------------------------------------------------------------
%   v       state vector (column vector of real numbers)
%   N       spatial resolution
%	symm    center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   f       the field in physical state (column vector of real numbers)


%%% REMARKS ---------------------------------------------------------------
%	1-  In forward integration of the KSE, the mean remains constant.
%       We assume the mean value is always zero, hence the mean value of
%       'f' is zero by construction.
%   2-  We do not fill the padded elements in Fourier space, hence no
%       additional dealiasing is needed.


function f = vector2field(v,N,symm)
    Nd = round(N/3);

    F = zeros(N,1);
    F(2:Nd+1) = v(1:Nd)*1j;
    
    if ~symm
        F(2:Nd+1) = F(2:Nd+1) + v(Nd+1:2*Nd);
    end

    F(N-Nd+1:N) = flip(conj(F(2:Nd+1)));
    
    f = ifft(F,'symmetric');
end
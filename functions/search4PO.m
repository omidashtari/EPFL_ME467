%%% DESCRIPTION -----------------------------------------------------------
%   shooting-based Newton iterations for computing periodic orbits (POs)


%%% INPUTS ----------------------------------------------------------------
%   u0      guessed periodic point (state vector)
%   T       guessed period of the PO
%   dt      reference time step for the required time marchings
%   L       domain length
%   N       spatial resolution
%   symm    center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   u_best  exact periodic point (state vector)
% 	T_best  exact period of the PO


function [u_best,T_best] = search4PO(u0,T,dt,L,N,symm)
    [~,k] = domain(L,N);
    U0 = fft(vector2field(u0,N,symm));
    
    V0 = -fft(ifft(U0,'symmetric').*ifft(complex(0,k).*U0,'symmetric'));
    V0 = dealiase(V0 + (k.^2 - k.^4).*U0);
    
    if symm
        V0 = complex(0,imag(V0));
    end
    
    v0 = field2vector(ifft(V0,'symmetric'),N,symm);
    
    options = optimoptions('fsolve','Display','iter','MaxIterations',75);
    [S,~,flag,~] = fsolve(@(X) recurrence(X),[u0;T],options);
    
    if flag ~= 1
        error('Search for unstable periodic orbit failed...');
    end
    
    u_best = S(1:end-1);
    T_best = S(end);
    
    function F = recurrence(X)
        u0_ = X(1:end-1);
        T_ = X(end);
        
        [uT_,~] = KSE_integrate(u0_,T_,dt,0,L,N,symm);
        F = [uT_-u0_; dot(u0_-u0,v0)];
    end
end
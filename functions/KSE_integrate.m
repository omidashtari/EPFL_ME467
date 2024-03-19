%%% DESCRIPTION -----------------------------------------------------------
%   time-marching the Kuramoto-Sivashinsky equation using ETDRK4 scheme
%   based on Kassam & Trefethen, SIAM J. Sci. Comput. 26, 1214-1233 (2005)


%%% INPUTS ----------------------------------------------------------------
%   u0          initial condition (column state vector)
%   T           integration time
%   dt_ref      reference time step size
%   dt_store    time intervals of storing (set to 0 to save 'u(T)' only)
%   L           domain length
%   N           spatial resolution
%   symm        center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   snapshots   stored snapshots (matrix with columns being state vectors)
%   t_vec       vector of instances associated with columns of 'snapshots'


%%% REMARKS ---------------------------------------------------------------
%   1-  Set 'dt_store = 0' to not save intermediate snapshots. In this case
%       'snapshots' contains only 'u(T)'.
%   2-  To recover the physical field corresponding to each column of the
%       'snapshots' matrix, use 'vector2field', and to transform a physical
%       field to the corresponding state vector use 'field2vector'.


function [snapshots,t_vec] = KSE_integrate(u0,T,dt_ref,dt_store,L,N,symm)
    %% assert data format
    if ~isreal(u0) || size(u0,2) ~= 1 || ~islogical(symm)
        error('ERROR: check type/size of inputs of KSE_integrate.')
    end
    
    %% adjust time step size
    N_steps = ceil(T/dt_ref);
    dt = T / N_steps;

    %% Fourier grid and the initial condition
    [~,k] = domain(L,N);
    v = fft(vector2field(u0,N,symm));

    %% pre-compute ETDRK4 scalars
    Linear = k.^2 - k.^4;
    
    E = exp(dt*Linear);
    E2 = exp(dt*Linear/2);
    
    M = 32;
    r = exp(1i*pi*((1:M)-0.5)/M);
    LR = dt*Linear(:,ones(M,1)) + r(ones(N,1),:);
    Q = dt*real(mean((exp(LR/2)-1)./LR,2));
    f1 = dt*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3,2));
    f2 = dt*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3,2));
    f3 = dt*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3,2));
    
    %% time-stepping loop:
    g = -0.5i*k;

    if dt_store == 0
        n = 1;
        Nt = 1;
    else
        n = floor(dt_store/dt);
        Nt = ceil(N_steps/n) + 1;
    end
    
    snapshots = zeros(length(u0),Nt);
    snapshots(:,1) = field2vector(ifft(v,'symmetric'),N,symm);
    
    t_vec = zeros(Nt,1);
    t_vec(1) = 0;

    t = 0;
    nt = 1;
    for q = 1:N_steps
        t = t+dt;
        
        Nv = g.*dealiase(fft(ifft(v,'symmetric').^2));
        
        a = E2.*v + Q.*Nv;
        Na = g.*dealiase(fft(ifft(a,'symmetric').^2));
        
        b = E2.*v + Q.*Na;
        Nb = g.*dealiase(fft(ifft(b,'symmetric').^2));
        
        c = E2.*a + Q.*(2*Nb-Nv);
        Nc = g.*dealiase(fft(ifft(c,'symmetric').^2));
        
        v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
        
        if symm
            v = complex(0,imag(v));
        end
        
        if dt_store ~= 0
            if rem(q,n) == 0
                nt = nt + 1;
                snapshots(:,nt) = field2vector(ifft(v,'symmetric'),N,symm);
                t_vec(nt) = t;
            end
        end
    end
    
    if t_vec(end) ~= T
        snapshots(:,end) = field2vector(ifft(v,'symmetric'),N,symm);
        t_vec(end) = T;
    end

    if dt_store == 0
        snapshots = field2vector(ifft(v,'symmetric'),N,symm);
        t_vec = T;
    end
end
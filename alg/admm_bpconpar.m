function [xsol, z, L1_v, L2_v, delta_v, sol_v, snr_v] = admm_bpconpar(y, epsilon, epsilons, A, At, T, W, Psi, Psit, param)
%
% [xsol, L1_v, L2_v, delta_v, sol_v, snr_v] = admm_bpconpar(y, epsilon, epsilons, A, At, T, W, Psi, Psit, param)
% implements Algorithm 1 as described in [1].
% The implementation simulates a distributed setup for parallel processing.
%
% Inputs:
% y{:} - the visibility data
% epsilon - global L2 bound
% epsilons - global L2 stop criterion bound
% A - function handle, linear operator modeling the measurement process,
%     FFT and zero padding 
% At - the adjoint of A
% T{:} - gridding matrix containing the interpolation kernels (and any other
%     modeled DDEs)
% W{:} - selection of the FFT coefficients associated with each block of
%        data
% Psi - function handle, prior regularisation function
% Psit - the adjoint of Psi
% param - configuration parameters
% 
% Outputs: 
% xsol - final solution
% L1_v - evolution of the L1 norm
% L2_v - evolution of the L2 norm
% delta_v - evolution of the relative solution variation
% sol_v - solutions at inerations probvided in param.sol_steps
% snr_v - evolution of the SNR
%
% Authors: Rafael Carrillo & Alexandru Onose
%
% [1] A. Onose, R. E. Carrillo, A. Repetti, J. D. McEwen, J.-P. Thiran,
% J.-C. Pesquet and Y. Wiaux - Scalable splitting algorithms for big-data
% interferometric imaging in the SKA era, MNRAS (2016) 
% doi:10.1093/mnras/stw1859 

% number of nodes
R = length(y);
P = length(Psit);

% oversampling vectorized data length
No = size(W{1}, 1);
[Noy, Nox] = size(A(At(zeros(No, 1))));

% number of pixels
N = numel(At(zeros(No, 1)));
[Ny, Nx] = size(At(zeros(No, 1)));

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'gamma'), param.gamma = 1e-2; end
if ~isfield(param, 'nu'), param.nu = 1; end

% Input arguments for prox L1
param_L1.Psi = Psi; 
param_L1.Psit = Psit; 
if isfield(param, 'nu_L1')
    param_L1.nu = param.nu_L1;
end
if isfield(param, 'tight_L1')
    param_L1.tight = param.tight_L1;
end
if isfield(param, 'max_iter_L1')
    param_L1.max_iter = param.max_iter_L1;
end
if isfield(param, 'rel_obj_L1')
    param_L1.rel_obj = param.rel_obj_L1;
else
    param_L1.rel_obj = param.rel_obj;
end
if isfield(param, 'weights')
    param_L1.weights = param.weights;
else
    param_L1.weights = 1;
end
if isfield(param, 'pos_L1')
    param_L1.pos = param.pos_L1;
else
    param_L1.pos = 0;
end
if isfield(param, 'verbose_L1')
    param_L1.verbose = param.verbose_L1;
else
    param_L1.verbose = param.verbose;
end

if ~isfield(param, 'sol_steps')
    param.sol_steps = inf;
else
    if param.sol_steps(end) ~= inf
        param.sol_steps = [param.sol_steps inf];
    end
end

% set up log variables
L1_v = zeros(param.max_iter, 1);
L2_v = zeros(param.max_iter, 1);

delta_v = zeros(param.max_iter, 1);
sol_steps = param.sol_steps;
sol_step_count = 1;
sol_v = zeros(length(sol_steps)-1, Ny, Nx);

snr_v = zeros(param.max_iter, 1);

%Initializations.

%Initial solution
if isfield(param,'initsol')
    xsol = param.initsol;
else
    xsol = zeros(size(At(zeros(size(W{1}))))); 
end

%Initial dual variables
if isfield(param, 'initz')
    z = param.initz;
else
    for k = 1:R
        z{k} =zeros(size(y{k}));
    end
end

%Initial residual and gradient (can be done in parallel)
res1 = zeros(R,1);
ns = A(xsol);
for k = 1:R
    %Residual
    res{k} = T{k}*ns(W{k}) - y{k};
    %Norm of residual
    res1(k) = norm(res{k});
    %Slack variable update
    s{k} = - z{k} - res{k};
    %Norm of slack variables
    normsl{k} = norm(s{k});
    %Projection onto the epsilon ball
    s{k} = min(epsilon(k)/normsl{k},1)*s{k};
    %Gradient formation
    nst = zeros(No, 1);
    nst(W{k}) = T{k}'*(z{k} + res{k} + s{k});
    g{k} = At(nst);
end

%Sum reduce to compute gradient
r = g{1};
for k = 2:R
    r = r + g{k};
end

%Flags initialization
dummy = Psit(xsol);
fval = sum(param_L1.weights(:).*abs(dummy(:))); 
flag = 0;

%Step sizes computation

%Step size primal 
mu = 1.0/param.nu;

%Step size for the dual variables
beta = 0.9;

epsilont = norm(epsilon);

%Main loop. Sequential.
for t = 1:param.max_iter
    
    tm = tic;
    
    %Gradient decend
    r = xsol - mu*r;
    
    prev_xsol = xsol;
    norm_prevsol = norm(prev_xsol(:));
    
    %Prox L1 norm (global solution)
    [xsol, fval] = solver_prox_L1(r, param.gamma*mu, param_L1);
    
    % solution relative change
    if (norm_prevsol == 0)
        rel_sol_norm_change = 1;
    else
        rel_sol_norm_change = norm(xsol(:) - prev_xsol(:))/norm_prevsol;
    end
    
    ns = A(xsol);
    %Parallel jobs
    for k = 1:R
        %Residual
        res{k} = T{k}*ns(W{k}) - y{k};
        %Norm of residual
        res1(k) = norm(res{k});
        %Slack variable update
        s{k} = - z{k} - res{k};
        %Norm of slack variables
        normsl{k} = norm(s{k});
        %Projection onto the epsilon ball
        s{k} = min(epsilon(k)/normsl{k},1)*s{k};
        %Lagrange multipliers update
        z{k} = z{k} + beta*(res{k} + s{k});
        %Gradient formation
        nst = zeros(No, 1);
        nst(W{k}) = T{k}'*(z{k} + res{k} + s{k});
        g{k} = At(nst);
    end
    
    tm = toc(tm);
    
    %Log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t);
        fprintf(' L1 norm = %e, rel solution norm change = %e\n', ...
            fval, rel_sol_norm_change);
        fprintf(' epsilon = %e, residual = %e\n\n', epsilont, norm(res1));
        for h = 1:R
            fprintf('  epsilon%i = %e, residual%i = %e\n\n', h, epsilon(h), h, res1(h));
        end
        
        fprintf('Time for iteration %i: %3.3f\n',t, tm);
        
    end
    
    if (param.verbose <= 0.5)
        fprintf('.\n');fprintf('\b');
        if mod(t, 50) == 0
            fprintf('\n');
        end
    end
    if (param.verbose >= 0.5)
        L1_v(t) = fval;
        L2_v(t) = norm(res1);
        
        delta_v(t) = rel_sol_norm_change;
        try 
            snr_v(t) = 20*log10(norm(param.im(:))/norm(param.im(:) - xsol(:)));
        catch
            snr_v(t) = 0;
        end
        if t == sol_steps(sol_step_count)
            sol_v(sol_step_count, :, :) = xsol;
            sol_step_count = sol_step_count + 1;
        end
    end
    %Global stopping criteria
    if (rel_sol_norm_change < param.rel_obj && prod(res1 <= epsilons))
        flag = 1;
        break;
    end
    
    %Sum reduce to compute gradient
    r = g{1};
    for k = 2:R
        r = r + g{k};
    end
          
end


%Final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf(' Objective function = %e\n', fval);
        fprintf(' Final residual = %e\n', res1);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf(' Objective function = %e\n', fval);
        fprintf(' Relative variation = %e\n', rel_sol_norm_change);
        fprintf(' Final residual = %e\n', res1);
        fprintf(' epsilon = %e\n', epsilon);
    end
end

% trim the log vectors to the size of the actual iterations performed
if (param.verbose >= 0.5)
    L1_v = L1_v(1:t);
    L2_v = L2_v(1:t);
    delta_v = delta_v(1:t);
    sol_v = sol_v(1:sol_step_count-1, :, :);
    snr_v = snr_v(1:t);
end

end


function [relerr_orig, relerr_sub, err_sub] = main(n, f, num_direction, num_prec, s, use_existing_respol)
%
% This is the main function of the first experiment.
%
% Input:
%
% - n: size of the problem (graph);
%
% - f: input vector for the computation. If empty, a random vector will be generated;
%
% - num_prec: flag indicating the number of preconditioners used to
% generate directions (reduced basis);
%       = 1 for 1 preconditioner. This is the common case;
%       = 2 for 2 preconditioners. Two preconditioners can improve accuracy;
%
% - s: the exponent in A^(-s). The default is set to be 0.5;
%
% - use_existing_respol: flag to decide whether to use provided residues and poles. If false, they will be computed using the OGA function.
%       = 1 for using provided residues and poles in the paper.
%       = 0 for using OGA function to compute residues and poles.
% Output:
%
% - relerr_orig = relative error of the original problem A^s u = f,
%                    where the exact solution is approximated by the rational approximation method. 
%
% - relerr_sub: relative error of each subproblem A_i u = f;
%
% - err_sub: absolute error of each subproblem A_i u = f.
%
% The function generates a graph and computes matrices based on the graph and input size 'n'. 
% It then performs preconditioning if required and computes the rational approximation for the function based on the residues and poles of the OGA rational approximant. 
% The relative and absolute errors are computed for each pole.



    % Load existing residues and poles or compute them using OGA
    if nargin < 6 || use_existing_respol
        [res, pol] = res_pols;
    else
        [res, pol, OGA_err] = OGA(s, 20);
    end

    % Generate graph and define matrices
    rng(100); % We use this seed in the experiment.
    L = generate_graph(n); 
    I = speye(n);
    A = L + (1 / n) * I;
    A = A / norm(A, inf); 

    % If input vector f is not provided, generate random vector
    if nargin < 2 || isempty(f)
        f = randn(n, 1);
    end

    nf = length(f); 
    npol = length(pol);
    x0 = sparse(nf, 1);
    B = A;

    % Perform preconditioning if required
    if nargin < 4 || (num_prec == 1)
        [~, P] = pcg(A + I, f, 1e-12, num_direction, B, x0);
    else
        pol_ = sort(pol);
        A_max = A + pol_(end) * I;
        [~, P1] = pcg(A + I, f, 1e-12, fix(num_direction / 2), B, x0);
        [~, P2] = pcg(A, f, 1e-12, fix(num_direction / 2), A_max, x0);
        P = [P1 P2];
    end

    % Precompute product matrices
    P = orth(P);
    PAPt = P' * A * P;      
    Pf = P' * f;  
    PAx0 = P' * A * x0;
    PMPt = P' * I * P;

    % Initialize error and result vectors
    err_sub = zeros(1, npol);
    relerr_sub = zeros(1, npol);
    u = zeros(n, 1);
    u_app = zeros(n, 1);
    
    % Set options for AMG
    option = struct('tol', 1e-10, 'printlevel', 0);

    % Loop over poles
    for k = 1:npol
        A_k = A + pol(k) * I;

        % Use AMG for large systems, otherwise solve directly
        if n > 2^10
            uk_exact = amg(A_k, f, option);
        else
            uk_exact = A_k \ f;
        end

        % Compute approximated solution
        y = (PAPt + pol(k) * PMPt) \ (Pf - PAx0);
        uk_app = P * y;

        % Compute solutions for A^s u = f. u is the rational approximation
        % to A^s = f. u_app is the approximation to u by RBMCG.
        u = u + res(k) * uk_exact; % u is the rational approximation solution to A^s u = f.
        u_app = u_app + res(k) * uk_app; % u_app is the approximation given by RCGBM. 

        % Compute errors for each subproblem A_i u = f
        err_sub(k) = norm(uk_exact - uk_app);    
        relerr_sub(k) = err_sub(k) / norm(uk_exact);
    end

    relerr_orig = norm(u - u_app) / norm(u);
end

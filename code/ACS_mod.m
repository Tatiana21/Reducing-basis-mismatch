clear all;

Q = 1;  % degree of overcompleteness
S = 6;  % sparsity
N = 256; 
M = 128; % number of measurements
dim = Q * N;
SNR = 40;% in dB
rng(9); %rng(5);
% A = randn(M,N); %sensing matrix

%generate testset
% [y, z, f_org, k_org, theta_org, phi_org] = generate_y(N,M,Q,S,SNR,A);
load demodata
%% ACS 
alpha = 0.1;
beta = 0.05;
TOL = 1e-4;
t = 1;
n_iter = 25; %temporarily fix number of iterations and calculate tolerance value
tol_val = zeros(1,n_iter);
rmse = zeros(1,n_iter);

theta = zeros(dim/2,1);
x_old = zeros(dim,1);

%function handles for computing psi and cost fn
psi_fn = @(theta) compute_psi(theta,Q,N);
f_fn = @(x, psi, lambda) norm_fn(x, psi, A, y, lambda);

lambda = alpha*(norm((A*psi_fn(theta))'*y,inf));
nonzeros = [];
%%
DEBIAS = 1;
xinit = (A*psi_fn(theta))'*y;
[x,x_debias] = GPSR_BB(y,A*psi_fn(theta),lambda,'Initialization',xinit,'Debias',DEBIAS,'Verbose',0);
if DEBIAS
    x_final = x;
    x       = x_debias;
end
xinit = x;

psi = psi_fn(theta);
fhat = psi*xinit;
%while (abs((f(x,psi_fn(theta)) - f(x_old,psi))) / f(x_old,psi)>TOL)
while(t <= n_iter)
    tic;
    tol_val(t)= (f_fn(x,psi_fn(theta), lambda) - f_fn(x_old,psi,lambda)) / f_fn(x_old,psi,lambda);
    rmse(t) = norm(f-psi*x)/norm(f);
    S_l = double([]);
    x_old = x;
    theta_old = theta;
        
    lambda = alpha*(norm((A*psi)'*y,inf));
    [x,x_debias] = GPSR_BB(y,A*psi,lambda,'Initialization',xinit,'Debias',DEBIAS,'Verbose',0); %update x
    if DEBIAS
        x_final = x;
        x       = x_debias;
    end
    th = 0.4;
    xinit = x;
    nonzeros = find(abs(x(2:end))>th)+1; % Call vector a nonzero if it is greater than a threshold - don't include DC (why I shift 2:end)
    nonzeros = [nonzeros;flipud(N-nonzeros+2)]; % Make sure we have cosine AND sine coefficients
    v = sort(nonzeros);                         % because if we perturb one we have to perturb the other!
    ind = find([1;diff(v)]&v);
    nonzeros = v(ind);
    toc;
    tic
%     theta = coordinate_GD(y, A, psi_fn, x, S_l, Q, N);
    psi = coordinate_GD_mod(y, A, x, nonzeros, N, psi);
    toc;
    t = t+1;
    fprintf(1,'\n',t);
end

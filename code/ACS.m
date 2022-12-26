clear all;

Q = 1;  % degree of overcompleteness
S = 6;  % sparsity
N = 256; 
M = 128; % number of measurements
dim = Q * N;
SNR = 40;% in dB
rng(9); %rng(5);
A = randn(M,N); %sensing matrix

%generate testset
[y, z, f, k_org, theta_org, phi_org] = generate_y(N,M,Q,S,SNR,A);
% load demodata
%% ACS 
alpha = 0.1;
beta = 0.05;
TOL = 1e-4;
t = 1;
n_iter = 10; %fix number of iterations and calculate tolerance value
tol_val = zeros(1,n_iter);
rmse = zeros(1,n_iter);

theta = zeros(dim/2,1);
x_old = zeros(dim,1);

%function handles for computing psi and cost fn
psi_fn = @(theta) compute_psi(theta,Q,N);
f_fn = @(x, psi, lambda) norm_fn(x, psi, A, y, lambda);

lambda = alpha*(norm((A*psi_fn(theta))'*y,inf));

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
%while (abs((f(x,psi_fn(theta)) - f(x_old,psi))) / f(x_old,psi)>TOL)
while(t <= n_iter)
    tic;
    tol_val(t)= abs((f_fn(x,psi_fn(theta), lambda) - f_fn(x_old,psi,lambda)) / f_fn(x_old,psi,lambda));
    rmse(t) = norm(z-psi*x)/norm(z);
    S_l = double([]);
    x_old = x;
    theta_old = theta;
    
    psi = psi_fn(theta);   
    lambda = alpha*(norm((A*psi)'*y,inf));
    [x,x_debias] = GPSR_BB(y,A*psi,lambda,'Initialization',xinit,'Debias',DEBIAS,'Verbose',0); %update x
    if DEBIAS
        x_final = x;
        x       = x_debias;
    end
    
    xinit = x;
    kappa = beta*(norm(x));
    for l=2:dim
       if (kappa<abs(x(l)) && l <= (dim / 2))
        S_l = [S_l l];
       elseif (kappa<abs(x(l)) && l >= (dim / 2))
       	S_l = [S_l (dim - l + 1)];
       end
    end
    
    v = sort(S_l);
    ind = find([1,diff(v)]&v);
    S_l = v(ind);       %support set 
    toc;
    tic
    theta = coordinate_GD(y, A, psi_fn, x, S_l, Q, N);
    toc;
    t = t+1;
    disp(t);
end
%% Plots
figure;plot(y);
hold on
plot(A*psi*x);
legend ('y','A*psi*x');
title('Plot of measured and recovered signal');

plvec = zeros(1,N);
plvec(k_org(1:3))= 0.5;
figure; plot(plvec(1:N))
hold on
plot(x)
axis([0, N, 0, inf]);
legend('original','recovered')
title('Frequency grid points')
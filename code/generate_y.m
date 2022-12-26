function [y,z,f, k, theta, phi] = generate_y(N,M,Q,S,SNR,A)
dim = Q*N;
k = randperm(floor(dim/2-1),S/2)+1;
theta = -1/(2*dim) + 1/(dim)*rand(dim/2,1);
phi = 2*pi*rand(S/2,1);   %phase
f = zeros(dim/2,1);         %frequencies
z = zeros(N,1);

for i=1:S/2
    f(i) = k(i)/(dim) + theta(k(i));
end
for i=1:N
    for s=1:S/2
        z(i) = sqrt(2/N)*cos(2*pi*f(s)*i + phi(s)) + z(i);
%         z(i) = sqrt(2/N)*cos(2*pi*(f(s))*i+phi(s)) + z(i);
    end
end

sigma = std(z)/(10^(SNR/20));
noise = sigma*randn(M,1);

y = A*z + noise;
end
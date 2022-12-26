function [ psi ] = compute_psi( theta, Q, N)
psi = zeros(N,Q*N);
for n=1:N
    for k=1:(Q*N/2)
        psi(n,k)= sqrt(2/N)*cos(2*pi*(n-1)*((k-1)/(Q*N) + theta(k)));
    end
    for k=(Q*N/2)+1:(Q*N)
        psi(n,k)= (-1)*sqrt(2/N)*sin(2*pi*(n-1)*((Q*N-k)/(Q*N) - theta(Q*N-k+1)));
    end
end
end


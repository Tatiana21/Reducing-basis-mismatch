function [psi,pertsol] = coordinate_GD_mod(y, A, x, supp, N, psi)
df=1/N;
options = optimset('TolX',df/10000,'TolCon',df/10000,'MaxIter',500,'Algorithm','interior-point','Display','off','GradObj','on','SubProblemAlgorithm','cg','Hessian','fin-diff-grads');

pertvec=zeros(N,1);
ub=df/2*ones(N,1);
lb=-df/2*ones(N,1);

offset = zeros(length(supp),N);
offset(length(supp)/2+1:end,:) = pi/2;

t1=repmat([0:1:N-1],length(supp),1);

for i = 1:length(supp)/2
    tmp=zeros(length(supp),1);
    tmp(i)=1;
    tmp(length(supp)-i+1) = -1;
    
    S=@(pvec)sqrt(2/N)*cos(2*pi*t1.*repmat(((supp-1)./(N)+pertvec(supp)+pvec*tmp)',N,1)'+offset)';
    Sprime=@(pvec2)(diag(tmp)*2*pi*sqrt(2/N)*(-1)*t1.*sin(2*pi*t1.*repmat(((supp-1)./(N)+pertvec(supp)+pvec2*tmp)',N,1)'+offset))';
    
    % solve for vector perturbations
    pertsol = fmincon(@(pvec)deal(  norm(y-A*S(pvec)*x(supp))^2, 2*(-A*Sprime(pvec)*x(supp))'*(y-A*S(pvec)*x(supp))),0.0,0,0,0,0,lb(supp(i)),ub(supp(i)),[],options); 
    
    lb(supp(i)) = lb(supp(i)) - pertsol; 
    ub(supp(i)) = ub(supp(i)) - pertsol;  % Update valid range for perturbations
        
    pertvec(supp(i)) = pertvec(supp(i)) + pertsol;  % Update the perturbation vector itself for the next go-round
    pertvec(N-supp(i)+2) = pertvec(N-supp(i)+2) - pertsol;
    
    psi(:,supp(i)) = sqrt(2/N)*cos(2*pi*[0:1:N-1]*((supp(i)-1)/(N)+pertvec(supp(i))));  % Update the basis vectors
    psi(:,N-supp(i)+2) = -sqrt(2/N)*sin(2*pi*[0:1:N-1]*((N-supp(i)+2-1)/(N)+pertvec(N-supp(i)+2)));
    
end
end
function [theta] = coordinate_GD(y, A, psi_fn, x, supp, Q, N)
% This function computes \theta that satisfies the relation
% \theta = arg min \norm{y - A * \Psi_theta * x}_2
% using coordinate gradient descent on the support set.
	
	dim = Q * N;
	theta = zeros(dim / 2,1);
	granularity = 1e-2 / dim; % quantization level for search
	scan_space = [- 1 / (2*dim) : granularity : 1 / (2*dim)];
    cost_l = zeros(size(scan_space));
    
	for l = supp
        l;
		% solve for \theta_l and update
		theta_temp = theta;
		for m = 1:size(scan_space,2)
			theta_temp(l) = scan_space(m);
			cost_l(m) = sum((y - (A * psi_fn(theta_temp) * x)).^2);
		end

		[~, ind] = min(cost_l);
		theta(l) = scan_space(ind);
	end

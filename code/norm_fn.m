function [ f_x_theta ] = norm_fn(x, psi, A, y, lambda )

f_x_theta = norm(y - A*psi*x)^2 + lambda*norm(x,1);

end


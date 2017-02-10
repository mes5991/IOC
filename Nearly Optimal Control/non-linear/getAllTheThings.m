function [mon, grad, hess, w, theta] = getAllTheThings(num_vars, num_monomials)

vars = sym('x', [1, num_vars]);
mon = monomials(vars, 1:num_monomials); %phi, or psi? 

grad = jacobian(mon, vars); %Gradient in Prof. Fu's paper.
hess = hessian(mon, vars); %How to do hessian of a vector?

w = sym('w', [1, num_vars]); %Weights for approximated value function, V
theta = sym('th', [1, num_vars]); %Weights for approximated cost function, l 



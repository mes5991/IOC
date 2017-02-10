function [ X ] = ode_fun( t,x)

X = [0;0];
u = u_nonlinear(x(1),x(2));
X(1) = x(1) + x(2) -x(1)*(x(1)^2 + x(2)^2);
X(2) = -x(1) + x(2) -x(2)*(x(1)^2 + x(2)^2) + u;

end

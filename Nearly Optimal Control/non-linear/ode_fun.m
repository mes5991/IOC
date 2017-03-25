function [ X ] = ode_fun( t,x,u_fun,noise_ratio)

X = [0;0];
u = u_fun(x(1),x(2));
X(1) = x(1) + x(2) -x(1)*(x(1)^2 + x(2)^2) +rand()*noise_ratio;
X(2) = -x(1) + x(2) -x(2)*(x(1)^2 + x(2)^2) + u +rand()*noise_ratio;

end


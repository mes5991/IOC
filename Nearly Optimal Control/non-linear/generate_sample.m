function [ X,u ] = generate_sample( start, tspan,u_func )
%GENERATE_SAMPLE Summary of this function goes here
%   Detailed explanation goes here

[t,X] = ode45(@(t,y)ode_fun(t,y,@(x1,x2)u_func(x1,x2)),[tspan],start);
u = zeros(size(t));
for i=1:length(t)
    u(i) = u_func(X(i,1),X(i,2));
end

end


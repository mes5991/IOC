function [ X,u ] = generate_sample( start, tspan,u_func,noise_ratio )
%GENERATE_SAMPLE Summary of this function goes here
%   Detailed explanation goes here

[t,X] = ode45(@(t,y)ode_fun(t,y,@(x1,x2)u_func(x1,x2),noise_ratio),[tspan],start);
u = zeros(size(tspan));
for i=1:length(tspan)
    u(i) = u_func(X(i,1),X(i,2));
end

end


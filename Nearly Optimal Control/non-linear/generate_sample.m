function [ X,u ] = generate_sample( start, tspan )
%GENERATE_SAMPLE Summary of this function goes here
%   Detailed explanation goes here
[t,X] = ode45(@(t,y)ode_fun(t,y),[tspan],start);
u = zeros(size(tspan));
for i=1:length(tspan)
    u(i) = u_nonlinear(X(i,1),X(i,2));
end

end


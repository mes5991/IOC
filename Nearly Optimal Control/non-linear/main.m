x = [0;1];

[t,X] = ode45(@(t,y)ode_fun(t,y),[0,30],x);
plot(t,X);
u = zeros(size(t));
figure();
for i=1:length(t)
    u(i) = u_nonlinear(X(i,1),X(i,2));
end
plot(t,u);
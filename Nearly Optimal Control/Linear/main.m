x = [1;1;1];

[t, X] = ode45(@(t,y)dynamics(t, y),[0,10],x);

x1 = X(:,1);
x2 = X(:,2);
x3 = X(:,3);
u = zeros(length(X),2);

for i = 1:length(X)
    u(i,1) = -3 * tanh((7.7*x1(i) + 2.44*x2(i) + 4.8*x3(i) ...
               +2.45*x1(i)^3 + 2.27*(x1(i)^2)*x2(i) + 3.7*x1(i)*x2(i)*x3(i) ...
               +0.71*x1(i)*x2(i)^2 + 5.8*(x1(i)^2)*x3(i) + 4.8*x1(i)*x3(i)^2 ...
               +0.08*x2(i)^3 + 0.6*(x2(i)^2)*x3(i) + 1.6*x2(i)*x3(i)^2 + 1.4*x3(i)^3)/3);
           
    u(i,2) = -20 * tanh((9.8*x1(i) + 2.94*x2(i) + 2.44*x3(i) ...
                -0.2*x1(i)^3 - 0.02*(x1(i)^2)*x2(i) + 1.42*x1(i)*x2(i)*x3(i) ...
                +0.12*x1(i)*x2(i)^2 + 2.3*(x1(i)^2)*x3(i) + 1.9*x1(i)*x3(i)^2 ...
                +0.02*x2(i)^3 + 0.23*(x2(i)^2)*x3(i) + 0.57*x2(i)*x3(i)^2 + 0.52*x3(i)^3)/20);
    
end


%%
figure(1)
plot(t, X)
grid on

figure(2)
plot(t,u)
grid on
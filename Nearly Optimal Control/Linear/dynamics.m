function [X] = dynamics(t, x)

x1 = x(1);
x2 = x(2);
x3 = x(3);

u1 = -3 * tanh((7.7*x1 + 2.44*x2 + 4.8*x3 ...
               +2.45*x1^3 + 2.27*(x1^2)*x2 + 3.7*x1*x2*x3 ...
               +0.71*x1*x2^2 + 5.8*(x1^2)*x3 + 4.8*x1*x3^2 ...
               +0.08*x2^3 + 0.6*(x2^2)*x3 + 1.6*x2*x3^2 + 1.4*x3^3)/3);
           
u2 = -20 * tanh((9.8*x1 + 2.94*x2 + 2.44*x3 ...
                -0.2*x1^3 - 0.02*(x1^2)*x2 + 1.42*x1*x2*x3 ...
                +0.12*x1*x2^2 + 2.3*(x1^2)*x3 + 1.9*x1*x3^2 ...
                +0.02*x2^3 + 0.23*(x2^2)*x3 + 0.57*x2*x3^2 + 0.52*x3^3)/20);
            
            
X = [2*x1 + x2 + x3;
     x1 - x2 + u2;
     x3 + u1];
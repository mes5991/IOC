function [ u ] = u_nonlinear( X1,X2)
%U_NONLINEAR Summary of this function goes here
%   Detailed explanation goes here
weights = [2.6,4.2,0.4,-4.0,-8.7,-8.9,-5.5,2.26,5.8,11,2.6,2.0,2.1,-0.5,-1.7,-2.71,-2.19,-.08,1.8,0.9];
x1 = X1;
x2= X2;
monomials = [x1; x2; x2^3; x1^3; x1^2*x2; x1*x2^2; x2^5; x1^5; x1^4*x2; ...
             x1^3*x2^2; x1^2*x2^3; x1*x2^4; x2^7; x1^7; x1^6*x2; x1^5*x2^2; ...
             x1^4*x2^3; x1^3*x2^4; x1^2*x2^5; x1*x2^6];
u = -tanh(weights*monomials);

u = subs(u);

end


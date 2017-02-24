%x = [-1:.1:1;-1:.1:1];
syms x1 x2
weights = [2.6,4.2,0.4,-4.0,-8.7,-8.9,-5.5,2.26,5.8,11,2.6,2.0,2.1,-0.5,-1.7,-2.71,-2.19,-.08,1.8,0.9];

mono = [x1; x2; x2^3; x1^3; x1^2*x2; x1*x2^2; x2^5; x1^5; x1^4*x2; ...
             x1^3*x2^2; x1^2*x2^3; x1*x2^4; x2^7; x1^7; x1^6*x2; x1^5*x2^2; ...
             x1^4*x2^3; x1^3*x2^4; x1^2*x2^5; x1*x2^6];
         
uF = matlabFunction(-tanh(weights*mono));

x = [0;1];

t = 0:0.01:30;


[X,u] = generate_sample([x(1);x(2)],t,@(x1,x2)uF(x1,x2));


%%
syms x1 x2
%h = [x1; x2; x2^3; x1^3; x1^2*x2; x1*x2^2; x2^5; x1^5; x1^4*x2; ...
%             x1^3*x2^2; x1^2*x2^3; x1*x2^4; x2^7; x1^7; x1^6*x2; x1^5*x2^2; ...
%             x1^4*x2^3; x1^3*x2^4; x1^2*x2^5; x1*x2^6];
variables = sym('a',[1,length(x)]);
h = monomials(variables,[1:7]);
hf = matlabFunction(h);
H = reshape(hf(X(:,1),X(:,2)),length(t),length(h))';


%%
lambda = 0.039
f = @(eta)abs(sum((tanh(eta'*H) - u).^2)) + lambda*eta'*eta
eta = lsqnonlin(f,zeros(size(h)))




control_law = @(x1,x2)tanh(eta'*hf(x1,x2));

t2 = [0:0.01:30];
[Xcomputed,ucomputed] = generate_sample([x(1);x(2)],t2,@(x1,x2)control_law(x1,x2));

close all
figure;
plot(t,u)
hold on
plot(t,tanh(eta'*H),'-.r')
plot(t2,ucomputed,'-.b')
hold off
legend('expert','fitted','computed')
figure
plot(t,X)
hold on
plot(t2,Xcomputed,'-.')
hold off
legend('x1 expert','x2 expert','x1 computed','x2 computed')
%x = [-1:.1:1;-1:.1:1];
syms x1 x2
weights = [2.6,4.2,0.4,-4.0,-8.7,-8.9,-5.5,2.26,5.8,11,2.6,2.0,2.1,-0.5,-1.7,-2.71,-2.19,-.08,1.8,0.9];

mono = [x1; x2; x2^3; x1^3; x1^2*x2; x1*x2^2; x2^5; x1^5; x1^4*x2; ...
             x1^3*x2^2; x1^2*x2^3; x1*x2^4; x2^7; x1^7; x1^6*x2; x1^5*x2^2; ...
             x1^4*x2^3; x1^3*x2^4; x1^2*x2^5; x1*x2^6];
         
uF = matlabFunction(-(weights*mono));

x = [0;1];

t = 0:0.01:69;


[X,u] = generate_sample([x(1);x(2)],t,@(x1,x2)uF(x1,x2));


%%
syms x1 x2
%h = [x1; x2; x2^3; x1^3; x1^2*x2; x1*x2^2; x2^5; x1^5; x1^4*x2; ...
%             x1^3*x2^2; x1^2*x2^3; x1*x2^4; x2^7; x1^7; x1^6*x2; x1^5*x2^2; ...
%             x1^4*x2^3; x1^3*x2^4; x1^2*x2^5; x1*x2^6];
variables = sym('a',[1,length(x)]);
h = monomials(variables,[1:8]);
hf = matlabFunction(h);
H = reshape(hf(X(:,1),X(:,2)),length(t),length(h))';


%%
lambda = 0.039
%lambda = 0.1
Dynamics = @(eta)abs(sum(((eta'*H) - u).^2)) + lambda*eta'*eta
eta = lsqnonlin(Dynamics,zeros(size(h)))





%%
control_law = @(x1,x2)(eta'*hf(x1,x2));
x = [-0,1]
t2 = [0:0.01:69];
[Xcomputed,ucomputed] = generate_sample([x(1);x(2)],t2,@(x1,x2)control_law(x1,x2));


%% Compute Lyapunov Control
syms x1 x2 real
vars = [x1;x2];
% Constructing the vector field dx/dt = f
Dynamics = [x1 + x2 - x1*(x1^2 + x2^2);
     -x1 + x2 - x2*(x1^2 + x2^2) + control_law(x1,x2)];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vars);

% =============================================
% The Lyapunov function V(x): 
mons = monomials([x1,x2],[1:8])
[prog,V] = sospolyvar(prog,mons,'wscoeff');

% =============================================
% Next, define SOSP constraints

% Constraint 1 : V(x) >= 0
prog = sosineq(prog,V);

% Constraint 2: (dV/dx)'*(f(x) + g(x)*u) <= 0
expr = -(gradient(V,[x1;x2])'*Dynamics+gradient(V,[x1;x2])'*Dynamics);
prog = sosineq(prog,expr);

% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Finally, get solution
V = sosgetsol(prog,V)


%% Formulate Control from CLF
% Assume sigma(x) = sqrt((Vx*(f(x)))^2 + Vx*g(x)'*R^(-1)*Vx*g(x))
Vx = gradient(V,[x1;x2]);
 f = [x1 + x2 - x1*(x1^2 + x2^2);
     -x1 + x2 - x2*(x1^2 + x2^2)];
 g = [0;1]
 LfV = Vx'*f;
 Lfg = Vx'*g;
 sigma = sqrt((LfV)^2 + (Lfg)^2);
phi = simplify(-(LfV + sigma)/(Lfg));
phi_fun = matlabFunction(phi)

%%
x = [0;1];
t2 = [0:0.001:69];
% XcomputedPhi = zeros([2,length(t2)]);
% ucomputedPhi = zeros([2,length(t2)]);
% ucomputedPhi(:,1) = phi_fun(XcomputedPhi(1,i),XcomputedPhi(2,i));
% for i=2:length(t2)
%     XcomputedPhi(:,i) = ode_fun(t,XcomputedPhi(:,i-1),@(x1,x2)phi_fun(x1,x2));
%     ucomputedPhi(:,i) = phi_fun(XcomputedPhi(1,i),XcomputedPhi(2,i));
% end
[XcomputedPhi,ucomputedPhi] = generate_sample([x(1);x(2)],t2,@(x1,x2)phi_fun(x1,x2));
%%
close all
figure;
plot(t,u)
hold on
plot(t,(eta'*H),'-.r')
plot(t2,ucomputed,'-.b')
plot(t2,ucomputedPhi,'--')
hold off
legend('expert','fitted','estimated','computed from Lyapunov')
figure
plot(t,X)
hold on
plot(t2,Xcomputed,'-.')
plot(t2,XcomputedPhi,'--')
legend('x1 expert','x2 expert','x1 estimated','x2 estimated','x1 computed from Lyapunov','x2 computed from Lyapunov')

%% compute sigma such that inf (Vx*(f(x)+g(x)*u) <= sigma(x)
sigma = -Vx'*f - Vx'*g;
l1 = (sigma - (Vx'*f))^2*(Vx'*g)^2
l1_fun = matlabFunction(l1)

l1Demonstrated = l1_fun(X(:,1),X(:,2))
l1Estimated = l1_fun(Xcomputed(:,1),Xcomputed(:,2))
l1ComputedPhi = l1_fun(XcomputedPhi(:,1),XcomputedPhi(:,2))
%%
figure
semilogy(t,l1Demonstrated)
hold on
semilogy(t2,l1Estimated,'-.')
semilogy(t2,l1ComputedPhi,'--')
grid on

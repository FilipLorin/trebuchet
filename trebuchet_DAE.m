close all; clc
% --------- %
% SOLVE DAE %
% --------- %
% https://www.mathworks.com/help/symbolic/solve-differential-algebraic-equations.html

%% SYMBOLIC PROBLEM DEFINITION

% Symbolic variables definition
syms th(t) phi(t) psi(t) lambda(t) M M14 M23 g I I1 I2 I14 I23 l2 l3

% Set of Differential Algebraic Equations (DAEs)
eqn1 = I*diff(th(t),2)+I14*sin(th(t)-phi(t))*(diff(phi(t),1))^2+I14*cos(th(t)-phi(t))*diff(phi(t),2)-I23*sin(th(t)-psi(t))*(diff(psi(t),1))^2-I23*cos(th(t)-psi(t))*diff(psi(t),2)+M*g*cos(th(t))-l2*cos(th(t))*lambda(t) == 0;
eqn2 = I1*diff(phi(t),2)-I14*sin(th(t)-phi(t))*(diff(th(t),1))^2+I14*cos(th(t)-phi(t))*diff(th(t),2)+M14*g*cos(th(t)) == 0;
eqn3 = I2*diff(psi(t),2)+I23*sin(th(t)-psi(t))*(diff(th(t),1)-diff(psi(t),1))*(diff(th(t),1))^2-I23*cos(th(t)-psi(t))*diff(th(t),2)+M23*g*cos(psi(t))+l3*cos(psi(t))*lambda(t) == 0;
eqn4 = -l2*cos(th(t))*diff(th(t),2)+l3*cos(psi(t))*diff(psi(t),2)+l2*(diff(th(t),1))^2*sin(th(t))-l3*(diff(psi(t),1))^2*sin(psi(t)) == 0;
eqns = [eqn1 eqn2 eqn3 eqn4];

% Variables
vars = [th(t); phi(t); psi(t); lambda(t)];
origVars = length(vars);

%% REDUCTION OF THE PROBLEM INDEX

% Incidence matrix (if one row is full of zeros there is a dof that is not
% necessary
Mat = incidenceMatrix(eqns,vars);

% Reduce order of the differential problem
[eqns,vars] = reduceDifferentialOrder(eqns,vars);

isLowIndexDAE(eqns,vars);

[DAEs,DAEvars] = reduceDAEIndex(eqns,vars);
[DAEs,DAEvars] = reduceRedundancies(DAEs,DAEvars);

isLowIndexDAE(DAEs,DAEvars);

pDAEs = symvar(DAEs);
pDAEvars = symvar(DAEvars);
extraParams = setdiff(pDAEs,pDAEvars);


%% SOLUTION OF THE PROBLEM
% Generation of the symbolic function handle

f = daeFunction(DAEs,DAEvars,I1, I2, I14, I23, I, M, M14, M23, g, l2, l3);

% Numerical parameters
m1 = 10;        % [kg]
m2 = 0.4;       % [kg]
ma = 2;         % [kg]
l1 = 0.6;       % [m]
l2 = 1.5;       % [m]
l3 = 0.7;       % [m]
l4 = 0.3;       % [m]
la = 0.5;       % [m]
h0 = 1.0;       % [m]
Ia = 2;         % [kgm^2]
g = 9.81;       % [ms^-2]

I1 = m1*l4*l4;
I2 = m2*l3*l3;
I14 = m1*l1*l4;
I23 = m2*l2*l3;
I = m1*l1*l1+m2*l2*l2+ma*la*la+Ia;
M = m1*l1-m2*l2-ma*la;
M14 = m1*l4;
M23 = m2*l3;

% Generation of numerical function handle
F = @(t,Y,YP) f(t,Y,YP,I1, I2, I14, I23, I, M, M14, M23, g, l2, l3);

% We generate the initial condition
th0 = asin(h0/l2);
phi0 = -pi/2;
psi0 = 0;
lambda0 = -M*g/l2;
Dth0 = 0;
Dphi0 = 0;
Dpsi0 = 0;

% Initial guess for the initial condition
y0est = [th0; phi0; psi0; lambda0; Dth0; Dphi0; Dpsi0];
yp0est = zeros(7,1);

% Find a consistent initial value set
opt = odeset('RelTol', 10.0^(-6),'AbsTol',10.0^(-7));
[y0,yp0] = decic(F,0,y0est,[],yp0est,[],opt);

% We solve the numerical problem
[tSol,ySol] = ode15i(F,[0 2],y0,yp0,opt);

% We extract the positions
thSol = ySol(:,1);
phiSol = ySol(:,2);
psiSol = ySol(:,3);

r1x = l1*cos(thSol)+l4*cos(phiSol);
r1y = l1*sin(thSol)+l4*sin(phiSol);
r2x = -l2*cos(thSol)+l3*cos(psiSol);
r2y = -l2*sin(thSol)+l3*sin(psiSol);
rax = -la*cos(thSol);
ray = -la*sin(thSol);

% Graphical representation
fid = 1;
figure(fid), fid=fid+1;
hold on
plot(tSol,ySol(:,1:origVars),'LineWidth',2)
for k = 1:origVars
  S{k} = char(DAEvars(k));
end
legend(S,'Location','Best')
grid on


figure(fid), fid=fid+1;
plot(r1x,r1y,'o');
hold on
plot(r2x,r2y,'s');
plot(rax,ray,'*');
hold off
grid on
axis equal
legend




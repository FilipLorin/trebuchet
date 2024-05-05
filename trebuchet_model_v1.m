%% Control variables
clc;

%%% Discrete time space
time_start = 0; %[s]
time_end = 2; %[s]
time_steps = 1000;

t_space = linspace(time_start,time_end,time_steps);
delta_t = (time_end-time_start)/time_steps;

%% Masses and Dimentions
m_1 = 10; %[kg] - counterweight mass
m_2 = 0.4; %[kg] - projectile mass
m_a = 2; %[kg] - arm mass

l_1 = 1; %[m] - arm length: central pivot to counterweight pivot
l_2 = 1.2; %[m] - arm length: central pivot to sling pivot
l_3 = 0.3; %[m] - sling length
l_4 = 0.5; %[m] - counterweight sling length
l_a = 0.2; %[m] - central pivot to arm COM (positive towards projectile)
h_0 = 0.7; %[m] - height: central pivot to floor
% NOTE - h_0 < l_2, else, please change initial values in solution

I_a = 2; %[kgm^2] - MOI of arm about COM

hook_angle = 0.7*pi; %[rad] - angle of release hook
prop_angle = -pi/2; %[rad] - absolute starting angle of counterweight arm

%% Initial optimisation by Krzysztof (overrides previous values)

l = 1.6; %fixed [m] - lenght of trebuchet is limited by the potential base contruciton dimensions
sigma = 0.0015; %fixed [m] - metal profile thickness
a = 0.050; %fixed [m] - main profile dimensions a,b
b = 0.030; %fixed [m]
ro = 7600; %fixed [kg/m^3] - density of our steel
m_support_1 = 0.23; %fixed [kg] - arm supports wieghts
m_support_2 = 0.2; %fixed [kg]
m_support_3 = 0.34; %fixed [kg]

l_1 = 0.8; %arm length: central pivot to counterweight pivot (0.1-0.8)[m] - we decided to place center of mass in front of pivot point so i'm not sure if we can move it to negatives
l_2 = l - l_1; % l_2 depends the lenghts of l_1 and l
l_3 = 0.55; %sling lenght (0.1-l_2)[m] - if sling lenght would be larger than l_2 our mechanism won't work anmore
l_4 = 0.05; %(0-0.3)[m] from the constructor side, its bad to make long counterweight mass arm
l_a = l_2 - l/2; %l_a is always on the middle point of our arm (assumption)
h_0 = 0.7; %(h_0 < l_2) [m] (as Filip sugested)

m_1 = 15; %fixed [kg] (Its going to be concrete I think)
m_2 = 0.045; %fixed (golf ball)[kg]
m_a = (2*l*a*sigma + 2*l*(b-2*sigma)*sigma)*ro + 2*m_support_1 + 2*m_support_2 + m_support_3; % (summ of all parts)

I_a = 1/12*m_a*l*l + 2*l*(b-2*sigma)*sigma*ro * (a-sigma)/2*(a-sigma)/2 + 2*0.3*0.3*m_support_1 + 2*0.14*0.14*m_support_2 + 0.2*0.2*m_support_3; % (based on simpe cad model) 


%% Preprocessing

%%% Constants
g = 9.81; %[m/s^2]

%%% Derived static values
M = m_1*l_1+m_2*l_2-m_a+l_a;
M_14 = m_1*l_4;
M_23 = m_2*l_3;

I = m_1*l_1*l_1+m_2*l_2*l_2+m_a*l_a*l_a+I_a;
I_1 = m_1*l_4*l_4;
I_2 = m_2*l_3*l_3;
I_14 = m_1*l_1*l_4;
I_23 = m_2*l_2*l_3;

%% Solution
%%% Euler method:
q = zeros(3, time_steps); % DOF angles: theta, phi, psi
q_dot = zeros(3, time_steps);
q_ddot = zeros(3, time_steps);

% Initial DOF values
q(:, 1) = [asin(h_0/l_2), prop_angle, 0];

release = time_steps; % indicator of frame on witch projectile was released
still_on_floor = 0;

for i = 2:time_steps
    M_matrix = [
        I,  I_14*cos(q(1, i-1)-q(3, i-1)),  -I_23*cos(q(1,i-1)-q(3, i-1));
        I_14*cos(q(1, i-1)-q(2, i-1)),  I_1,    0;
        -I_23*cos(q(1, i-1)-q(3, i-1)), 0, I_2
        ];
    g_vec = [
        -I_14*sin(q(1, i-1)-q(2, i-1))*q_dot(2, i-1)^2+I_23*sin(q(1,i-1)-q(3, i-1))*q_dot(3, i-1)^2-M*g*cos(q(1, i-1));
        I_14*sin(q(1, i-1)-q(2, i-1))*q_dot(1, i-1)^2-M_14*g*cos(q(2, i-1));
        -I_23*sin(q(1, i-1)-q(3, i-1))*q_dot(1, i-1)^2-M_23*g*cos(q(3,i-1))
        ];
    
    if h_0-l_2*sin(q(1, i-1))+l_3*sin(q(3, i-1)) > 0 % projectile off the floor
        if q(3, i-1) + 2*pi < q(1, i-1)+hook_angle % projectile released
            release = i-1;
            break
        else % unconstrained movement
            q_ddot(:, i) = M_matrix\g_vec;
        end
    else % projectile on floor
        still_on_floor = i-1;
        partial_fq = [-l_2*cos(q(1,i-1)), 0, l_3*cos(q(3,i-1))]';
        A = [ M_matrix, partial_fq;
            partial_fq', 0 ];
        s = -l_2*q_dot(1, i-1)^2*sin(q(1, i-1))+l_3*q_dot(3, i-1)^2*sin(q(3, i-1));
        b = [g_vec; s];
        x = A\b;
        q_ddot(:,i) = x(1:3);
    end
    q_dot(:,i) = q_dot(:,i-1)+q_ddot(:,i-1)*delta_t;
    q(:,i) = q(:,i-1) + q_dot(:,i-1)*delta_t + q_ddot(:,i-1)*delta_t^2/2;
end

%%% Inital results
disp("Projectile left ground at [s] t=")
disp(t_space(still_on_floor))
disp("Projectile launched at [s] t=")
disp(t_space(release))
r_2_dot = [l_2*sin(q(1, release))*q_dot(1, release)-l_3*sin(q(3, release))*q_dot(3, release);
            -l_2*cos(q(1, release))*q_dot(1, release)+l_3*cos(q(3, release))*q_dot(3, release)];
disp("Release velocity [m/s]:")
disp(r_2_dot)

%% Postprocessing
close all;

%%% Mass position vectors
r_1 = [l_1*cos(q(1,1:release))+l_4*cos(q(2,1:release)); h_0+l_1*sin(q(1,1:release))+l_4*sin(q(2,1:release))];
r_2 = [-l_2*cos(q(1,1:release))+l_3*cos(q(3,1:release)); h_0-l_2*sin(q(1,1:release))+l_3*sin(q(3,1:release))];
r_a = [-l_a*cos(q(1,1:release)); h_0-l_a*sin(q(1,1:release))];


figure % Mass trajectory plot
plot(r_1(1,:), r_1(2,:), "--x")
hold on
grid on
axis equal
plot(r_2(1,:), r_2(2,:), "--x")
plot(r_a(1,:), r_a(2,:), "--x")
line([-2, 2], [0, 0])
line([0, 0], [0, h_0])
legend("r_1", "r_2", "r_a", location='southeast')


% figure % plot of DOF angles in time
% plot(t_space(1:release), 180/pi*q(1,1:release))
% hold on
% plot(t_space(1:release), 180/pi*q(2, 1:release))
% plot(t_space(1:release), 180/pi*q(3, 1:release))
% grid on
% ticks = -240:30:60;
% yticks(ticks)
% legend("$\Theta$", "$\varphi$", "$\psi$", interpreter="latex", location="southwest")

% figure % Projectile y/time plot
% plot(t_space(1:release), r_2(2,:))
% grid on
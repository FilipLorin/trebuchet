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

% The main values to optimalize are l_1, l_3 (sling lenght), l_4 (counterweight arm), h_0.
% Also we need to optimize hook angle to get right release (I suggest release angle = pi/4);

% l_4 makes the smallest changes to velocity vector, so we can ignore it if we want to.

% For me our goal should be to maxmize absolute velocity with release angle around pi/4 (my personal
% minimum is 28m/s (100km/h)).

% On default (given above) values with hook_angle = 0.72*pi; we got [22.35 22.7]' output projectile velocity
% which gives us absolute velocity around 31.8 m/s.
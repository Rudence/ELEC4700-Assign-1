% Creator: Rudi Hidvary 
% Student Number: 101037815
% Class: ELEC 4700 
% Document: Assignment 1

clear
clc

% Constants for Model 
m0 = 9.11e-31; % electron mass (kg)
k = 1.381e-23; % boltzmans constant 

% Model Parameters
length = 200e-9;        % size of simulation in x dierection (m)
height = 100e-9;        % size of simulation in y direction (m)
temperature = 300;      % temperature in kelvin
me = 0.26*m0;           % Effective mass of an electorn in our simulation
e_num = 400;             % Number of electrons in the simulation 
simlength = 2000;        % Sests the number of iterations the simulation undergoes
graph_pause = 1;
sim_pause = 0.00001;


% Initial Calculations 
% Question 1.a THERMAL VELOCITY
thermal_velocity = sqrt((2*k*temperature)/me) % velocity in (m/s)

% Initializing the Simulation Parameters 
initial_xposition = length*rand(e_num,1); % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
initial_yposition = height*rand(e_num,1); % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation
theta = 2*pi*rand(e_num,1); % Initializes a vector of random angles the size of the number of electrons
initial_xvelocity = thermal_velocity.*cos(theta).*ones(e_num,1); % Sets the initial velocity as the thermal temperature  
initial_yvelocity = thermal_velocity.*sin(theta).*ones(e_num,1); % 

figure(1)
plot(initial_xposition, initial_yposition, 'o')
title('Initial Particle Positions')
xlabel('X Position (m)')
ylabel('Y Position (m)')
axis([0 200e-9 0 100e-9])
grid on
pause(graph_pause)

timestep = 1e-15; % Timestep is the amount of time between each interval of the calculations 

new_xposition = initial_xposition; % Sets the 
new_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

temp = [];

% Simulation loop that continually updates the simulation at each timestep
% and calculates the 
for time = 1:simlength 
        
new_xposition = new_xposition + new_xvelocity*timestep;
new_yposition = new_yposition + new_yvelocity*timestep;

% Boundary Conditions being imposed 
overboundx = new_xposition > 200e-9;
underboundx = new_xposition < 0;
overboundy = new_yposition > 100e-9;
underboundy = new_yposition < 0;
new_xposition(overboundx) = new_xposition(overboundx) - 200e-9;  
new_xposition(underboundx) = new_xposition(underboundx) + 200e-9;
new_yvelocity(overboundy) = -new_yvelocity(overboundy);
new_yvelocity(underboundy) = -new_yvelocity(underboundy);

% Boundary Conditions for Bottleneck region


% Plotting the updating positions 
% Question 1.c.i 2D PLOT OF TRAJECTORIES
figure(2)
plot(new_xposition,new_yposition,'ro')
title('Simulation')
xlabel('Distance (nm)')
ylabel('Distance (nm)')
grid on
axis([0 200e-9 0 100e-9]) 
pause(sim_pause)

averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
temp(time) = (averageVel*me)/(2*k);

end
hold off

time = 1:simlength;

% Question 1.c.ii TEMPERATURE PLOT
figure(3) 
plot(time,temp)
title('')
xlabel('')
ylabel('')
grid on
axis([time(1) time(end) 0 max(temp)+100])
pause(graph_pause)

% Question 1.b MEAN FREE PATH
% Mean Free Path Calculation 
Tmn = 0.2e-12; % Mean time between collisions 
MFP = averageVel*Tmn % Mean didtance travelled before collision occurs








% Question 2.a
% Initial Calculations 
% Question 1.a THERMAL VELOCITY
thermal_velocity = sqrt((2*k*temperature)/me) % velocity in (m/s)
std_thermal_velocity = 0.1*thermal_velocity;
random_velocity = normrnd(thermal_velocity,std_thermal_velocity,[e_num,1]);
bin_num = 5;
figure(4)
velocity_hist = histogram(random_velocity,bin_num);
title('Thermal Velocity Distribution')
xlabel('Random Thermal Velocity (m/s)')
ylabel('Number of Particles Within Range')
pause(graph_pause)

% Initializing the Simulation Parameters 
initial_xposition = length*rand(e_num,1); % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
initial_yposition = height*rand(e_num,1); % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation
theta = 2*pi*rand(e_num,1); % Initializes a vector of random angles the size of the number of electrons
initial_xvelocity = random_velocity.*cos(theta).*ones(e_num,1); % Sets the initial velocity as a value from the distribution that gives thermal velocity
initial_yvelocity = random_velocity.*sin(theta).*ones(e_num,1); % 

figure(5)
plot(initial_xposition, initial_yposition, 'o')
title('Initial Particle Positions')
xlabel('X Position (m)')
ylabel('Y Position (m)')
axis([0 200e-9 0 100e-9])
grid on
pause(graph_pause)

timestep = 1e-15; % Timestep is the amount of time between each interval of the calculations 

new_xposition = initial_xposition; % Sets the 
new_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

temp = [];

% Simulation loop that continually updates the simulation at each timestep
% and calculates the 
for time = 1:simlength 

% Electron scattering and reevaluation of velocity
rand_threshold = rand(e_num,1);
Pscatter = (1-exp(-(timestep/Tmn)))
    for index = 1:e_num
        if rand_threshold(index) < Pscatter 
            theta = 2*pi*rand(1);
            new_velocity = normrnd(thermal_velocity,std_thermal_velocity,[1,1]);
            new_xvelocity(index) = cos(theta)*new_velocity;  
            new_yvelocity(index) = sin(theta)*new_velocity;
        end
    end
    
new_xposition = new_xposition + new_xvelocity*timestep;
new_yposition = new_yposition + new_yvelocity*timestep;

% Boundary Conditions being imposed 
overboundx = new_xposition > 200e-9;
underboundx = new_xposition < 0;
overboundy = new_yposition > 100e-9;
underboundy = new_yposition < 0;
new_xposition(overboundx) = new_xposition(overboundx) - 200e-9;  
new_xposition(underboundx) = new_xposition(underboundx) + 200e-9;
new_yvelocity(overboundy) = -new_yvelocity(overboundy);
new_yvelocity(underboundy) = -new_yvelocity(underboundy);

% Plotting the updating positions 
% Question 1.c.i 2D PLOT OF TRAJECTORIES
figure(6)
plot(new_xposition,new_yposition,'ro')
title('Simulation')
xlabel('Distance (nm)')
ylabel('Distance (nm)')
grid on
axis([0 200e-9 0 100e-9]) 
pause(sim_pause)

averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
temp(time) = (averageVel*me)/(2*k);



end
hold off

time = 1:simlength;

% Question 1.c.ii TEMPERATURE PLOT
figure(7) 
plot(time,temp)
title('')
xlabel('')
ylabel('')
grid on
axis([time(1) time(end) 0 max(temp)+100])
pause(graph_pause)

% Question 1.b MEAN FREE PATH
% Mean Free Path Calculation 
Tmn = 0.2e-12; % Mean time between collisions 
MFP = averageVel*Tmn % Mean didtance travelled before collision occurs

















% % Question 3.a 2D PLOT OF TRAJECTORIES 
% 
% % Bottleneck Characteristics
% % These are all centred around the middle, the width is the size of the
% % neck in x position and height fixes the hegith of the neck
% bottle_width = 40e-9; % Remember that the max original width is 200-e9 m
% bottle_height = 40e-9; % Remember the max original height is 100e-9 m
% 
% % Initializing the Simulation Parameters 
% % No electrons start in the block regions 
% initial_xposition = (length-bottle_width)*rand(e_num,1); % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
% initial_yposition = (height-bottle_height)*rand(e_num,1); % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation
% 
% move_right = (initial_xposition > 0.5*(length-bottle_width)) & ((initial_yposition > 0.5*(height+bottle_height)) | (initial_yposition < 0.5*(height-bottle_height)))
% initial_xposition(move_right) = initial_xposition(move_right) + bottle_width;
% move_up = (initial_yposition > 0.5*(height-bottle_height))
% initial_yposition(move_up) = initial_yposition(move_up) + bottle_height;
% theta = 2*pi*rand(e_num,1); % Initializes a vector of random angles the size of the number of electrons
% initial_xvelocity = thermal_velocity.*cos(theta).*ones(e_num,1); % Sets the initial velocity as the thermal temperature  
% initial_yvelocity = thermal_velocity.*sin(theta).*ones(e_num,1); % 
% 
% figure(1)
% plot(initial_xposition, initial_yposition, 'o')
% title('Initial article Positions')
% xlabel('X Position (m)')
% ylabel('Y Position (m)')
% axis([0 200e-9 0 100e-9])
% grid on
% 
% timestep = 1e-15; % Timestep is the amount of time between each interval of the calculations 
% 
% new_xposition = initial_xposition; % Sets the 
% new_yposition = initial_yposition;
% new_xvelocity = initial_xvelocity;
% new_yvelocity = initial_yvelocity;
% 
% temp = [];
% 
% % Simulation loop that continually updates the simulation at each timestep
% % and calculates the 
% for time = 1:simlength 
%         
% new_xposition = new_xposition + new_xvelocity*timestep;
% new_yposition = new_yposition + new_yvelocity*timestep;
% 
% % Boundary Conditions being imposed 
% overboundx = new_xposition > 200e-9;
% underboundx = new_xposition < 0;
% overboundy = new_yposition > 100e-9;
% underboundy = new_yposition < 0;
% new_xposition(overboundx) = new_xposition(overboundx) - 200e-9;  
% new_xposition(underboundx) = new_xposition(underboundx) + 200e-9;
% new_yvelocity(overboundy) = -new_yvelocity(overboundy);
% new_yvelocity(underboundy) = -new_yvelocity(underboundy);
% 
% % Boundary Conditions for Bottleneck region
% 
% 
% % Plotting the updating positions 
% % Question 1.c.i
% figure(2)
% plot(new_xposition,new_yposition,'ro')
% title('Simulation')
% xlabel('Distance (nm)')
% ylabel('Distance (nm)')
% grid on
% axis([0 200e-9 0 100e-9]) 
% pause(0.01)
% 
% averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
% temp(time) = (averageVel*me)/(2*k);
% 
% end
% hold off
% 
% time = 1:simlength;
% 
% % Question 1.c.ii
% figure(3) 
% plot(time,temp)
% title('')
% xlabel('')
% ylabel('')
% grid on
% axis([time(1) time(end) 0 max(temp)+10])





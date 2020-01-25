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
e_num = 50;             % Number of electrons in the simulation 
simlength = 100;        % Sests the number of iterations the simulation undergoes

% Initial Calculations 
thermal_velocity = sqrt((k*temperature)/me) % velocity in (m/s)

% Initializing the Simulation Parameters 
initial_xposition = length*rand(e_num,1) % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
initial_yposition = height*rand(e_num,1) % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation
theta = 2*pi*rand(e_num,1) % Initializes a vector of random angles the size of the number of electrons
initial_xvelocity = thermal_velocity.*cos(theta).*ones(e_num,1) % Sets the initial velocity as the thermal temperature  
initial_yvelocity = thermal_velocity.*sin(theta).*ones(e_num,1) % 

figure(1)
plot(initial_xposition, initial_yposition, 'o')
title('Initial article Positions')
xlabel('X Position (m)')
ylabel('Y Position (m)')
axis([0 200e-9 0 100e-9])
grid on

timestep = 1e-15; % Timestep is the amount of time between each interval of the calculations 

new_xposition = initial_xposition; % Sets the 
new_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

temp = []

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

% Plotting the updating positions 
figure(2)
plot(new_xposition,new_yposition,'ro')
title('Simulation')
xlabel('Distance (nm)')
ylabel('Distance (nm)')
grid on
axis([0 200e-9 0 100e-9]) 
pause(0.01)

averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
temp(time) = (averageVel*me)/(2*k);

end
hold off

time = 1:simlength

figure(3) 
plot(time,temp)
title('')
xlabel('')
ylabel('')
grid on
axis([time(1) time(end) 0 max(temp)+10])

% Mean Free Path Calculation
Tmn = 0.2e-12; % Mean time between collisions 
MFP = averageVel*Tmn % Mean didtance travelled before collision occurs





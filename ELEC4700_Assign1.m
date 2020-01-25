% Creator: Rudi Hidvary 
% Student Number: 101037815
% Class: ELEC 4700 
% Document: Assignment 1


clear
clc

m0 = 9.11e-31; % electron mass (kg)
me = 0.26*m0;
length = 200e-9; % size of box (m)
height = 100e-9; % size of box (m)
temperature = 300; % temperature in kelvin
k = 1.381e-23; % boltzmans constant 
thermal_velocity = sqrt((k*temperature)/me) % velocity in (m/s)

e_num = 8


initial_xposition = 200*rand(e_num,1)
initial_yposition = 100*rand(e_num,1)
theta = 2*pi*rand(e_num,1)
initial_xvelocity = thermal_velocity.*cos(theta).*ones(e_num,1)
initial_yvelocity = thermal_velocity.*sin(theta).*ones(e_num,1)

figure(1)
plot(initial_xposition, initial_yposition, 'o')
axis([0 200 0 100])
grid on

timestep = 1/100000;

new_xposition = initial_xposition;
new_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

temp = []

for time = 1:20
        
new_xposition = new_xposition + new_xvelocity*timestep;
new_yposition = new_yposition + new_yvelocity*timestep;

% Boundary Conditions being imposed 
overboundx = new_xposition > 200;
underboundx = new_xposition < 0;
overboundy = new_yposition > 100;
underboundy = new_yposition < 0;
new_xposition(overboundx) = new_xposition(overboundx) - 200;  
new_xposition(underboundx) = new_xposition(underboundx) + 200;
new_yvelocity(overboundy) = -new_yvelocity(overboundy);
new_yvelocity(underboundy) = -new_yvelocity(underboundy);

% Plotting the updating positions 
figure(2)
plot(new_xposition,new_yposition,'ko')
title('Simulation')
xlabel('Distance (nm)')
ylabel('Distance (nm)')
grid on
axis([0 200 0 100]) 
pause(0.1)


averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
temp(time) = (averageVel*me)/(2*k);

figure(3) 
hold on
plot(time,temp)
grid on
% axis([0 time 0 200])



end
hold off

% Mean Free Path Calculation
Tmn = 0.2e-12; % Mean time between collisions 
MFP = averageVel*Tmn % Mean didtance travelled before collision occurs





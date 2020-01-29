% Creator: Rudi Hidvary 
% Student Number: 101037815
% Class: ELEC 4700 
% Document: Assignment 1

% QUESTIONS
% Q1: how do i get the trajectories not just the current values 
% Q2: mean free path for part 2 and how to get the time between collisions 
% Q3: How do i fix time on the temp graphs 

clear
clc

% Constants for Model 
m0 = 9.11e-31; % electron mass (kg)
k = 1.381e-23; % boltzmans constant 

%Model Parameters
length = 200e-9;        % size of simulation in x dierection (m)
height = 100e-9;        % size of simulation in y direction (m)
temperature = 300;      % temperature in kelvin
me = 0.26*m0;           % Effective mass of an electorn in our simulation
e_num = 500;             % Number of electrons in the simulation 
simlength = 200;        % Sests the number of iterations the simulation undergoes
graph_pause = 0.1;
sim_pause = 0.0000001;


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

initial_velocity = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
initial_temp = (initial_velocity*me)/(2*k);
temp = [initial_temp];

% Simulation loop that continually updates the simulation at each timestep
% and calculates the 
for time = 1:simlength 
        
averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2)); % Calculates the average velocity of the particles at given timestep
temp(time) = (averageVel*me)/(2*k); % Using the avaerge velocity, find the temperature of the system at this giventimestep

    
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

end
hold off

time = 1:simlength;

% Question 1.c.ii TEMPERATURE PLOT
figure(3) 
plot(time*timestep,temp)
title('Simulation Temperature Over Time')
xlabel('Time (s)')
ylabel('Simulation Temperature')
grid on
pause(graph_pause)

% Question 1.b MEAN FREE PATH
% Mean Free Path Calculation 
Tmn = 0.2e-12; % Mean time between collisions 
MFP_Q1 = averageVel*Tmn % Mean didtance travelled before collision occurs using the given mean time












%------------------------------------------------------------------------------------------------------------------------------------------------
% Question 2.a

% Initial Calculations 
% The thermal velocity needs to be made random and assigned to each
% variable. The random distribution is plotted in a histogram which shows
% the randomly generated distribution
thermal_velocity = sqrt((2*k*temperature)/me); % velocity in (m/s)
distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
random_velocity = maxwell_boltzmann_dist;
bin_num = 15;
figure(8)
velocity_hist = histogram(random_velocity,bin_num);
title('Thermal Velocity Distribution')
xlabel('Random Thermal Velocity (m/s)')
ylabel('Number of Particles Within Range')
grid on
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

initial_velocity = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
initial_temp = (initial_velocity*me)/(2*k);
temp = [initial_temp];

% Scattering Equations
Pscatter = (1-exp(-(timestep/Tmn))); % Used to compare to random variable for scattering chance

% Mean Free Path and Mean Collision Time Equations
MFP = [];
TMN = [];

% Simulation loop that continually updates the simulation at each timestep
% and calculates the 
for time = 1:simlength 

% Electron scattering and reevaluation of velocity
rand_threshold = rand(e_num,1); % Sets a vector of random numbers for each electron
distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
new_velocity = maxwell_boltzmann_dist;
    for index = 1:e_num
        if rand_threshold(index) < Pscatter 
            theta = 2*pi*rand(1);
            new_xvelocity(index) = cos(theta)*new_velocity(e_num);  
            new_yvelocity(index) = sin(theta)*new_velocity(e_num);
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
MFP(time) = averageVel*Tmn; % Uses the given Tmn to find the MFP for this simulation
TMN(time) = MFP_Q1/averageVel; % Uses the previously calculated MFP in part 1 to find the average time of collisions

end
hold off

% Question 1.c.ii TEMPERATURE PLOT
time = 1:simlength;

figure(7) % Plotting the Temperature of the system
plot(time*timestep,temp)
title('Average Temperature Over Time')
xlabel('Time (s)')
ylabel('Simulation Temperature (K)')
grid on
pause(graph_pause)

% Question 1.b MEAN FREE PATH
% Mean Free Path Calculation  
Tmn = 0.2e-12; % Mean time between collisions 
MFP_average = mean(MFP) % Mean distance travelled before collision occurs
TMN_average = mean(TMN) % Mean time between collisions calculated using a known 













%------------------------------------------------------------------------------------------------------------------------------------------------
% Question 3.a

% Initial Calculations 
% The thermal velocity needs to be made random and assigned to each
% variable. The random distribution is plotted in a histogram which shows
% the randomly generated distribution.
thermal_velocity = sqrt((2*k*temperature)/me)       ; % velocity in (m/s)
distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
random_velocity = maxwell_boltzmann_dist;
bin_num = 15;
figure(8)
velocity_hist = histogram(random_velocity,bin_num);
title('Thermal Velocity Distribution')
xlabel('Random Thermal Velocity (m/s)')
ylabel('Number of Particles Within Range')
grid on
pause(graph_pause)

% Setting whether the boundaries are specular = 0 or diffusive = 1
boundary_type = 0; % Specular
% boundary_type = 1; % Diffusive


% Initializing the Simulation Parameters 
initial_xposition = length*rand(e_num,1); % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
initial_yposition = height*rand(e_num,1); % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation

% Reinitializing electrons that spawned in restricted regions
% Length of the channel or bottleneck and height are determined here
Lbottle = 80e-9;
Hbottle = 20e-9;
% Checking the boundaries and reinitializing the ones that are restricted
inboundx = (initial_xposition < ((length/2)+(Lbottle/2))) & (initial_xposition > ((length/2)-(Lbottle/2)));
inboundy = (initial_yposition < ((height/2)-(Hbottle/2))) | (initial_yposition > ((height/2)+(Hbottle/2)));
inbound = inboundx & inboundy;
while(max(inbound) > 0)
initial_xposition(inbound) = rand(size(initial_xposition(inbound),1),1)*length;
initial_yposition(inbound) = rand(size(initial_yposition(inbound),1),1)*height;
%initial_xposition(inbound) = rand(size(inbound,1),1)*length;
%initial_yposition(inbound) = rand(size(inbound,1),1)*height;
inboundx = (initial_xposition < ((length/2)+(Lbottle/2))) & (initial_xposition > ((length/2)-(Lbottle/2)));
inboundy = (initial_yposition < ((height/2)-(Hbottle/2))) | (initial_yposition > ((height/2)+(Hbottle/2)));
inbound = inboundx & inboundy;
end


theta = 2*pi*rand(e_num,1); % Initializes a vector of random angles the size of the number of electrons
initial_xvelocity = random_velocity.*cos(theta).*ones(e_num,1); % Sets the initial velocity as a value from the distribution that gives thermal velocity
initial_yvelocity = random_velocity.*sin(theta).*ones(e_num,1); % 

figure(9)
plot(initial_xposition, initial_yposition, 'o')
title('Initial Particle Positions')
xlabel('X Position (m)')
ylabel('Y Position (m)')
axis([0 200e-9 0 100e-9])
grid on
pause(graph_pause)

timestep = 1e-15; % Timestep is the amount of time between each interval of the calculations 

old_xposition = initial_xposition; % Sets the 
old_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

initial_velocity = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
initial_temp = (initial_velocity*me)/(2*k);
temp = [initial_temp];

% Scattering Equations
Pscatter = (1-exp(-(timestep/Tmn))); % Used to compare to random variable for scattering chance

% Mean Free Path and Mean Collision Time Equations
MFP = [];
TMN = [];

% Simulation loop that continually updates the simulation at each timestep
% and calculates the 
for time = 1:simlength 

% Electron scattering and reevaluation of velocity
rand_threshold = rand(e_num,1); % Sets a vector of random numbers for each electron
distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
new_velocity = maxwell_boltzmann_dist;
    for index = 1:e_num
        if rand_threshold(index) < Pscatter 
            theta = 2*pi*rand(1);
            new_xvelocity(index) = cos(theta)*new_velocity(index);  
            new_yvelocity(index) = sin(theta)*new_velocity(index);
        end
    end
    
new_xposition = old_xposition + new_xvelocity*timestep;
new_yposition = old_yposition + new_yvelocity*timestep;

% Boundary Conditions being imposed 
overboundx = new_xposition > 200e-9;
underboundx = new_xposition < 0;
overboundy = new_yposition > 100e-9;
underboundy = new_yposition < 0;
new_xposition(overboundx) = new_xposition(overboundx) - 200e-9;  
new_xposition(underboundx) = new_xposition(underboundx) + 200e-9;
new_yvelocity(overboundy) = -new_yvelocity(overboundy);
new_yvelocity(underboundy) = -new_yvelocity(underboundy);

% Restricted Region Boundary Cases
overhorizontal =  (new_xposition > ((length/2)-(Lbottle/2))) & (new_xposition < ((length/2)+(Lbottle/2))) & ((new_yposition < ((height/2)-(Hbottle/2))) | (new_yposition > ((height/2)+(Hbottle/2))));
previous_left = (old_xposition < ((length/2)-(Lbottle/2)));
previous_right = (old_xposition > ((length/2)+(Lbottle/2)));
previous_in = (old_xposition > ((length/2)-(Lbottle/2))) & (old_xposition < ((length/2)+(Lbottle/2)));
% if particles come from left and go over to restricted region, flip the
% velocities 
% Diffusive Boundary (Random Generated new velocity)
if (boundary_type == 1)
    distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
    distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
    maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
    new_velocity = maxwell_boltzmann_dist;
    
    theta = 2*pi*rand(e_num,1);
    % For the electron to the left of the boundary
    new_xvelocity(overhorizontal & previous_left) = -abs(cos(theta(overhorizontal & previous_left))*new_xvelocity(overhorizontal & previous_left));
    new_yvelocity(overhorizontal & previous_left) = sin(theta(overhorizontal & previous_left))*new_yvelocity(overhorizontal & previous_left);
    % for the electrons to the right of the boundary
    new_xvelocity(overhorizontal & previous_right) = abs(cos(theta(overhorizontal & previous_left))*new_xvelocity(overhorizontal & previous_left));
    new_yvelocity(overhorizontal & previous_right) = sin(theta(overhorizontal & previous_left))*new_yvelocity(overhorizontal & previous_left);
    % for theelectron sin the tunnel region
    new_xvelocity(overhorizontal & previous_in) = cos(theta(overhorizontal & previous_left))*new_xvelocity(overhorizontal & previous_left);
    new_yvelocity(overhorizontal & previous_in) = sin(theta(overhorizontal & previous_left))*new_yvelocity(overhorizontal & previous_left);
    
end
else 
    new_xvelocity(overhorizontal & (previous_left | previous_right)) = -new_xvelocity(overhorizontal & (previous_left | previous_right));
    new_yvelocity(overhorizontal & previous_in) = -new_yvelocity(overhorizontal & previous_in);
end

% Plotting the updating positions 
% Question 3.a 2D PLOT OF TRAJECTORIES

% ATTEMPT AT TRAJECTORIES 
% plotx = [old_xposition new_xposition]
% plotx = plotx'
% ploty = [old_yposition new_yposition]
% ploty = ploty'
figure(10)
plot(plotx,ploty,'ro')
title('Simulation')
xlabel('Distance (nm)')
ylabel('Distance (nm)')
grid on
axis([0 200e-9 0 100e-9]) 
pause(sim_pause)

averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
temp(time) = (averageVel*me)/(2*k);
MFP(time) = averageVel*Tmn; % Uses the given Tmn to find the MFP for this simulation
TMN(time) = MFP_Q1/averageVel; % Uses the previously calculated MFP in part 1 to find the average time of collisions

old_xposition = new_xposition;
old_yposition = new_yposition;

end
hold off

% Question 1.c.ii TEMPERATURE PLOT
time = 1:simlength;

figure(11) % Plotting the Temperature of the system
plot(time*timestep,temp)
title('Average Temperature Over Time')
xlabel('Time (s)')
ylabel('Simulation Temperature (K)')
grid on
pause(graph_pause)

% Question 1.b MEAN FREE PATH
% Mean Free Path Calculation  
Tmn = 0.2e-12; % Mean time between collisions 
MFP_average = mean(MFP) % Mean distance travelled before collision occurs
TMN_average = mean(TMN) % Mean time between collisions calculated using a known 



figure(12)
hist3([new_xposition new_yposition],[10,10])
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');





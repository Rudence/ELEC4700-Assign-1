% Creator: Rudi Hidvary 
% Student Number: 101037815
% Class: ELEC 4700 
% Document: Assignment 1

% QUESTIONS
% Q1: how do i get the trajectories not just the current values 
% Q2: mean free path for part 2 and how to get the time between collisions 
% Q3: How do i fix time on the temp graphs 

% Clears all of the variables, the command line, and closes all previous
% graphs 
close all
clear
clc

% Constants for Model 
m0 = 9.11e-31; % electron mass (kg)
k = 1.381e-23; % boltzmans constant 

% Model Parameters
% ------------------------------------------------------------------------
% These parameters can be changed to see the effect they have on the
% simulation 
length = 200e-9;        % size of simulation in x dierection (m)
height = 100e-9;        % size of simulation in y direction (m)
temperature = 300;      % temperature in kelvin
me = 0.26*m0;           % Effective mass of an electron in our simulation
e_num = 10000;          % Number of electrons in the simulation 
simlength = 1000;       % Sets the number of iterations the simulation undergoes, 1 interation is 1 femtosecond
graph_pause = 0;        % Length graph is presented in a figure 
sim_pause = 0;          % Sets the length of time that the simulation pauses for at each step 
% the bin number calculation gives a reasonable bin number for any amount
% of particles over 30, 
bin_num =  10; %ceil(3*log10(e_num));           % Histogram bin number 
% The bin number for the amount of bins in the electron density plot for
% part 3
bin_num_3D = 20; 
% Choose the type of boundary conditions for the third simulation 
boundary_type =  0; % 0 is for specular 1 for diffusive boudary type
%-------------------------------------------------------------------------

% Initial Calculations 
% Question 1.a THERMAL VELOCITY
thermal_velocity = sqrt((2*k*temperature)/me); % velocity in (m/s)
fprintf('QUESTION 1\n')
fprintf('The thermal velocity of all of the particles Vth = %d',thermal_velocity);
fprintf('(m/s) \n')

% Initializing the Simulation Parameters 
initial_xposition = length*rand(e_num,1); % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
initial_yposition = height*rand(e_num,1); % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation
theta = 2*pi*rand(e_num,1); % Initializes a vector of random angles the size of the number of electrons
initial_xvelocity = thermal_velocity.*cos(theta).*ones(e_num,1); % Sets the initial velocity as the thermal temperature  
initial_yvelocity = thermal_velocity.*sin(theta).*ones(e_num,1); % 

% Plotting the initial positions of the electrons as they are distributed
% over the surface
figure(1)
plot(initial_xposition, initial_yposition, 'ko')
title('Initial Particle Positions')
xlabel('X Position (m)')
ylabel('Y Position (m)')
axis([0 200e-9 0 100e-9])
grid on
pause(graph_pause) % use a pause to see the graph as the simulatoin is running 

% Timestep is the amount of time between each interval of the calculations 
timestep = 1e-15; % 1 femtosecond 

new_xposition = initial_xposition; % Sets the 
new_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

initial_velocity = (mean(initial_xvelocity.^2)) + (mean(initial_yvelocity.^2));
initial_temp = (initial_velocity*me)/(2*k);
temp = [initial_temp];

% Simulation loop that continually updates the simulation at each timestep
% and calculates the new positions of the particles by using time and
% velocity. At each iteration the average velocity is found and 
for time = 1:simlength 
    % The average velocity and the temperature of the system at each time
    % step are calculated to be used in plotting and finding the mean free
    % path given the average time between collisions 
    averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2)); 
    temp(time) = (averageVel*me)/(2*k); 
    
    % Updating the new positions using the velocity
    new_xposition = new_xposition + new_xvelocity*timestep;
    new_yposition = new_yposition + new_yvelocity*timestep;

    % Boundary Conditions being imposed for particle wrap around and
    % reflectoin at the tep and bottom of the simulation surface
    % Logical index parameters that check if a praticle has left the bounds
    % and how these cases are handled
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
    figure(2)
    plot(new_xposition,new_yposition,'ro')
    title('Simulation')
    xlabel('Distance (nm)')
    ylabel('Distance (nm)')
    grid on
    axis([0 200e-9 0 100e-9]) 
    pause(sim_pause)
end
time = 0:simlength;
temp = [initial_temp temp];
% Question 1.c.ii TEMPERATURE PLOT
% Plotting the temperature as a function of time using the timestep for the
% x axis to get the actual time in the simulation. 
figure(3) 
plot(time*timestep,temp)
title('Simulation Temperature Over Time')
xlabel('Time (s)')
ylabel('Simulation Temperature (Kelvin)')
grid on
pause(graph_pause)

% Question 1.b MEAN FREE PATH
% Mean Free Path Calculation 
Tmn = 0.2e-12; % Mean time between collisions 
MFP_Q1 = averageVel*Tmn; % Mean didtance travelled before collision occurs using the given mean time
fprintf('The average distance a particle travels before a collision occurs MFP =  %d',MFP_Q1)
fprintf('(m) \n')
%-------------------------------------------------------------------------------------------------------------------------------------------



%------------------------------------------------------------------------------------------------------------------------------------------------
% Question 2.a

% Initial Calculations 
% The thermal velocity needs to be made random and assigned to each
% variable. The random distribution is plotted in a histogram which shows
% the randomly generated distribution. This is achieved by generating two 
% random distibutions and finding the square of the sum of squares to get 
% the maxwell-boltzmann distribution 
thermal_velocity = sqrt((2*k*temperature)/me); % velocity in (m/s)
distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
random_velocity = maxwell_boltzmann_dist;

% Plot of the histogram of the distributed velocities 
figure(4)
histogram(random_velocity,bin_num);
title('Thermal Velocity Distribution')
xlabel('Random Thermal Velocity (m/s)')
ylabel('Number of Particles Within Range')
grid on
pause(graph_pause)

% Initializing the Simulation Parameters 
% Setting the initial x and y velocities to values from the randomly
% generated distribution
initial_xposition = length*rand(e_num,1); % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
initial_yposition = height*rand(e_num,1); % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation
theta = 2*pi*rand(e_num,1); % Initializes a vector of random angles the size of the number of electrons
initial_xvelocity = random_velocity.*cos(theta).*ones(e_num,1); % Sets the initial velocity as a value from the distribution that gives thermal velocity
initial_yvelocity = random_velocity.*sin(theta).*ones(e_num,1); % 

% Initial position plot 
figure(5)
plot(initial_xposition, initial_yposition, 'ko')
title('Initial Particle Positions')
xlabel('X Position (m)')
ylabel('Y Position (m)')
axis([0 200e-9 0 100e-9])
grid on
pause(graph_pause)

timestep = 1e-15; % Timestep is the amount of time between each interval of the calculations 

new_xposition = initial_xposition;
new_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

initial_velocity = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
initial_temp = (initial_velocity*me)/(2*k);
temp = [initial_temp];

% Scattering Equation to get a probability value to be compared to 
Pscatter = (1-exp(-(timestep/Tmn))); % Used to compare to random variable for scattering chance

% Mean Free Path and Mean Collision Time Equation collision number
% initializatoin
collision_num = 0; 

% Question 2.b
% Simulation loop that continually updates the simulation at each timestep
% and keeps track of the number of collisoins that occur due to scattering
% and records the temperature of the system at each timestep.
for time = 1:simlength 
    % Electron scattering and reevaluation of velocity distribution at each
    % time step to get more random values than using the same distribution
    % that is outside of the loop
    rand_threshold = rand(e_num,1); % Sets a vector of random numbers for each electron to be compared to the scattering value
    distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
    distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
    maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
    new_velocity = maxwell_boltzmann_dist;
    
    % using a for loop to check each of the electrons in the system to see
    % if scattering occured by using a randomly generated number for each
    % of the particles, then rethermalizing the velocites if scatrering has occured 
    for index = 1:e_num
        if rand_threshold(index) < Pscatter 
            theta = 2*pi*rand(1);
            new_xvelocity(index) = cos(theta)*new_velocity(index);  
            new_yvelocity(index) = sin(theta)*new_velocity(index);
            %Counts the number of collisions that have occured to be used
            %to calculate the mean time between collisions 
            collision_num = collision_num + 1; 
        end
    end
    
    % Updating the new positions at each timestep
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

    % Plotting the updating positions of the particle trajectories
    % Question 2.b 2D PLOT OF TRAJECTORIES
    figure(6)
    plot(new_xposition,new_yposition,'ro')
    title('Simulation')
    xlabel('Distance (nm)')
    ylabel('Distance (nm)')
    grid on
    axis([0 200e-9 0 100e-9]) 
    pause(sim_pause)
    
    % finding the average velocities of all of the particles
    averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
    % Collect value of temperature over the simulation length
    temp(time) = (averageVel*me)/(2*k);
end

% Question 2.c TEMPERATURE PLOT
time = 0:simlength;
temp = [initial_temp temp];
% Plotting the Temperature of the system
figure(7)
plot(time*timestep,temp)
title('Average Temperature Over Time')
xlabel('Time (s)')
ylabel('Simulation Temperature (K)')
grid on
pause(graph_pause)

% Calculating the average amount of time between collisions by using the
% total time the simulation ran for and the number of particles in the
% simulation divided by the amount of collisions to get the average time
% betwen collisions
calculated_Tmn_Q2 = (timestep*simlength*e_num)/(collision_num);
% Once the average time between colliosns is calculated, then the mean free
% path can be found by mulitpling by the average velocity of all of the
% particles in the simulation
calculated_MFP_Q2 = calculated_Tmn_Q2*averageVel;

% Calculating the percent error of the mean time and mean free path
percent_error_Tmn = (100*abs(calculated_Tmn_Q2-Tmn))/Tmn;
percent_error_MFP = (100*abs(calculated_MFP_Q2-MFP_Q1))/MFP_Q1;

% Printing the answer statments 
fprintf('\nQUESTION 2')
fprintf('\nThe average time between collisions was found to be %d',calculated_Tmn_Q2)
fprintf(' (s)\nComparing this to the given value of %d (s)',Tmn)
fprintf('\nThe percent error of the average time is %f percent',percent_error_Tmn)
fprintf('\nThe mean free path was calculated to be %d',calculated_MFP_Q2)
fprintf(' (m)\nComparing this to the previously calculated value of %d (m)\n',MFP_Q1)
fprintf('The percent error of the MFP is %f percent\n',percent_error_MFP)
%-------------------------------------------------------------------------------------------------------------------------------------------------



%------------------------------------------------------------------------------------------------------------------------------------------------
% Question 3.a

% Initial Calculations 
% The thermal velocity needs to be made random and assigned to each
% variable. The random distribution is plotted in a histogram which shows
% the randomly generated distribution.
thermal_velocity = sqrt((2*k*temperature)/me); % velocity in (m/s)
distribution1 = randn(e_num,1)*(thermal_velocity/sqrt(2));
distribution2 = randn(e_num,1)*(thermal_velocity/sqrt(2));
maxwell_boltzmann_dist = sqrt((distribution1.^2)+(distribution2.^2));
random_velocity = maxwell_boltzmann_dist;

% Plotting the histogram of the randomly generated distribution 
figure(8)
histogram(random_velocity,bin_num);
title('Thermal Velocity Distribution')
xlabel('Random Thermal Velocity (m/s)')
ylabel('Number of Particles Within Range')
grid on
pause(graph_pause)

% Setting whether the boundaries are specular = 0 or diffusive = 1
% boundary_type = 0; % Specular
% boundary_type = 1; % Diffusive

% Initializing the Simulation Parameters 
initial_xposition = length*rand(e_num,1); % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
initial_yposition = height*rand(e_num,1); % Sets the initial y position as a vector of randomly selected numbers over the length of the simulation

% Reinitializing electrons that spawned in restricted regions
% Length of the channel or bottleneck and height are determined here
Lbottle = 80e-9;
Hbottle = 20e-9;

% Checking the boundaries and marking the particles that have spawned
% inside the restrcted region.
inboundx = (initial_xposition < ((length/2)+(Lbottle/2))) & (initial_xposition > ((length/2)-(Lbottle/2)));
inboundy = (initial_yposition < ((height/2)-(Hbottle/2))) | (initial_yposition > ((height/2)+(Hbottle/2)));
inbound = inboundx & inboundy;

% Using a while loop to reinitialize any of the particles that have spawned  in the restriced region 
while(max(inbound) > 0)
    initial_xposition(inbound) = rand(size(initial_xposition(inbound),1),1)*length;
    initial_yposition(inbound) = rand(size(initial_yposition(inbound),1),1)*height;
    % recheck the boundaries to see if there are still particles that have
    % spawned inside the restriced regions. 
    inboundx = (initial_xposition < ((length/2)+(Lbottle/2))) & (initial_xposition > ((length/2)-(Lbottle/2)));
    inboundy = (initial_yposition < ((height/2)-(Hbottle/2))) | (initial_yposition > ((height/2)+(Hbottle/2)));
    inbound = inboundx & inboundy;
end

% Initialize random velocities for each of the particles given the random
% distribution that was generated
theta = 2*pi*rand(e_num,1); % Initializes a vector of random angles the size of the number of electrons
initial_xvelocity = random_velocity.*cos(theta).*ones(e_num,1); % Sets the initial velocity as a value from the distribution that gives thermal velocity
initial_yvelocity = random_velocity.*sin(theta).*ones(e_num,1); % 

% Plots the initial positions of the particles 
figure(9)
plot(initial_xposition, initial_yposition, 'ko')
title('Initial Particle Positions')
xlabel('X Position (m)')
ylabel('Y Position (m)')
axis([0 200e-9 0 100e-9])
grid on
pause(graph_pause)

timestep = 1e-15; % Timestep is the amount of time between each interval of the calculations 

old_xposition = initial_xposition;  
old_yposition = initial_yposition;
new_xvelocity = initial_xvelocity;
new_yvelocity = initial_yvelocity;

initial_velocity = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
initial_temp = (initial_velocity*me)/(2*k);
temp = [initial_temp];

% Scattering Equations
Pscatter = (1-exp(-(timestep/Tmn))); % Used to compare to random variable for scattering chance
collision_num = 0;

% Simulation loop that continually updates the simulation at each timestep
% and calculates the average velocity
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
            %Counts the number of collisions that have occured to be used
            %to calculate the mean time between collisions 
            collision_num = collision_num + 1;
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
        new_xvelocity(overhorizontal & previous_left) = -abs(cos(theta(overhorizontal & previous_left)).*new_velocity(overhorizontal & previous_left));
        new_yvelocity(overhorizontal & previous_left) = sin(theta(overhorizontal & previous_left)).*new_velocity(overhorizontal & previous_left);
        % for the electrons to the right of the boundary
        new_xvelocity(overhorizontal & previous_right) = abs(cos(theta(overhorizontal & previous_right)).*new_velocity(overhorizontal & previous_right));
        new_yvelocity(overhorizontal & previous_right) = sin(theta(overhorizontal & previous_right)).*new_velocity(overhorizontal & previous_right);
        % for the electrons in the tunnel region
        new_xvelocity(overhorizontal & previous_in) = cos(theta(overhorizontal & previous_in)).*new_velocity(overhorizontal & previous_in);
        new_yvelocity(overhorizontal & previous_in) = sin(theta(overhorizontal & previous_in)).*new_velocity(overhorizontal & previous_in);
        % To boot any stray particles that may be stuck
        new_xposition(overhorizontal & ~previous_left & ~previous_right & ~previous_in) = 
    else 
        new_xvelocity(overhorizontal & (previous_left | previous_right)) = -new_xvelocity(overhorizontal & (previous_left | previous_right));
        new_yvelocity(overhorizontal & previous_in) = -new_yvelocity(overhorizontal & previous_in);
    end

    % Plotting the updating positions 
    % Question 3.a 2D PLOT OF TRAJECTORIES
    figure(10)
    plot(new_xposition,new_yposition,'ro')
    title('Simulation')
    xlabel('Distance (nm)')
    ylabel('Distance (nm)')
    grid on
    axis([0 200e-9 0 100e-9]) 
    pause(sim_pause)

    averageVel = (mean(new_xvelocity.^2)) + (mean(new_yvelocity.^2));
    temp(time) = (averageVel*me)/(2*k);

    old_xposition = new_xposition;
    old_yposition = new_yposition;

end

% Question 1.c.ii TEMPERATURE PLOT
time = 0:simlength;
temp = [initial_temp temp];

figure(11) % Plotting the Temperature of the system
plot(time*timestep,temp)
title('Average Temperature Over Time')
xlabel('Time (s)')
ylabel('Simulation Temperature (K)')
grid on
pause(graph_pause)

% Calculating the average amount of time between collisions by using the
% total time the simulation ran for and the number of particles in the
% simulation divided by the amount of collisions to get the average time
% betwen collisions
calculated_Tmn_Q3 = (timestep*simlength*e_num)/(collision_num);
% Once the average time between colliosns is calculated, then the mean free
% path can be found by mulitpling by the average velocity of all of the
% particles in the simulation
calculated_MFP_Q3 = calculated_Tmn_Q3*averageVel;

% Calculating the percent error fom the correct answer
percent_error_Tmn = (100*abs(calculated_Tmn_Q3-Tmn))/Tmn;
percent_error_MFP = (100*abs(calculated_MFP_Q3-MFP_Q1))/MFP_Q1;

% Printing the answer statments 
fprintf('\nQUESTION 3')
fprintf('\nThe average time between collisions was found to be %d',calculated_Tmn_Q3)
fprintf(' (s)\nComparing this to the given value of %d (s)',Tmn)
fprintf('\nThe percent error of the average time is %f percent',percent_error_Tmn)
fprintf('\nThe mean free path was calculated to be %d',calculated_MFP_Q3)
fprintf(' (m)\nComparing this to the preeviously calculated value of %d (m)\n',MFP_Q1)
fprintf('The percent error of the MFP is %f percent\n',percent_error_MFP)

% Plotting a 3D histogram for the elecron density over the surface of the simulation    
figure(12)
hist3([new_xposition new_yposition],[bin_num_3D bin_num_3D])
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('Electron Density Visual Map')
xlabel('Distance in X Direction (m)')
ylabel('Distance in Y Direction (m)')
zlabel('number of Electron in Area')
pause(graph_pause)

% Calculates the velocity of each of the particles so thaty they can be
% made into a temperature map.
particleVel = sqrt((new_xvelocity.^2) + (new_yvelocity.^2));
particleTemp = (particleVel.*me)./(2*k);

% Plots the colour spectrum for the particles for their calculated
% temperature and a bar to help dicern the temperature a=of a given
% particle
figure(13)
scatter(new_xposition,new_yposition,12,particleTemp)
colorbar
grid on
title('Coloured Temperature Map of Particles')
xlabel('Distance in X Direction (m)')
ylabel('Distance in Y Direction (m)')
axis([0 200e-9 0 100e-9])
pause(graph_pause)


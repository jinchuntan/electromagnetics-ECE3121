% Written by Tan Jin Chun
% Last Modified: 4/9/2023
% Code for the Simulation of Water Fall

clear all;close all;clc

% Plotting the Electric Field
% Entering a zero charge or omitting a charge will result in charges with 
% higher index numbers being neglected. (Input Charges)

% Go to Efieldplot.m file
Efieldplot();

% Plotting the equipotential and the 3D Surface Plot
global Vinitdisp X Y;
[Vinitdisp]=setVBCs;

[V,it,error]=laplacesolv();
V = flipud(V);

% Define constants (For the Simulation)
epsilon_0 = 8.854e-12;   % Vacuum permittivity (F/m)
rho_water = 1000;       % Density of water (kg/m^3)
r_droplet = 1e-4;       % Radius of the droplet (m) - just an example size
V1 = 120;                % Applied voltage across electrodes (V)
V2 = 0;
V3 = 120;
V4 = 120;
d = 1e-2;               % Distance between electrodes (m)
alpha_water = 1.35e-40; % Polarizability of water (F.m^2)
mass_water = 18.01528 / (6.022140857*10^23); % mass of water molecule divided by the Avogadro's Number

% Introducing Gravity
g = 9.81;
initial_height = 0;

% Constants for the Coulomb's law calculation
k = 8.9875e9; % Coulomb's constant (N.m^2/C^2)

% 1. Model the electric field between two parallel plates
E = V / d;  % Electric field (V/m) - for parallel plates

% 2. Compute the force on a water droplet
%    We will consider the dielectrophoretic force
F_de = 2 * pi * epsilon_0 * alpha_water * r_droplet^3 * E^2;


% Example charges on each electrode to create an electric field
% Area of the capacitance plate is assumed to be 
A = 0.02*0.02; % 0.02m * 0.02m, distance between the plate is assumed to be 0.01m
C1 = epsilon_0*A/0.005;
C2 = epsilon_0*A/0.005;
C3 = epsilon_0*A/0.005;
C4 = epsilon_0*A/0.005;
Q1 = C1 * V1;
Q2 = C2 * V2;
Q3 = C3 * V3;
Q4 = C4 * V4;
Q = [Q1, -Q2, Q3, -Q4];  % Charges on the four electrodes (C)

% Variables for the theoretical force acting on the water molecules
dE_dx = 10^3; % Assume to be approximately 10^3
alpha = 10^-24; % The polarizability of water
mu = 6.2*10^-30; % The dipole moment of a water molecule 
T = 368; % Assuming that the water temperatuere is near boiling point

% Positions of the electrodes (just example placements)
position_value = 0.02;
electrode_positions = [position_value, position_value; -position_value, position_value; -position_value, -position_value; position_value, -position_value];  % In meters

% Function to compute electric field at position 'pos' due to a point charge
computeE = @(q, pos, droplet_pos) k * q / norm(pos - droplet_pos)^3 * (pos - droplet_pos);

% 3. Simulate the motion of the droplet using Newton's 2nd law (F = ma)
%    For simplicity, we'll do a very basic time-stepping simulation
dt = 1e-3;   % Time step (s)
t_end = 1;   % Total simulation time (s)

% positions = zeros(t_end/dt, 1);  % Store positions over time
positions = zeros(t_end/dt, 2);  % Store positions over time in 2D
velocities = zeros(t_end/dt, 2); % Store velocities over time

% Initial conditions
positions(1) = 0;   % Initial position
velocities(1) = 0;  % Initial velocity

% Time-stepping loop
for t = 2:length(positions)

    % E Total
    E_total = [0, 0];  % Initialize total electric field to zero
    
    % Superposition: Add up electric field contributions from each electrode
    droplet_pos = positions(t-1, :);
    for i = 1:4
        E_from_electrode = computeE(Q(i), electrode_positions(i, :), droplet_pos);
        E_total = E_total + E_from_electrode;
    end

    % Compute DEP force based on E_total
    F_de = 2 * pi * epsilon_0 * alpha_water * r_droplet^3 * norm(E_total)^2;
    % F_de = theoreticalForce(E_total, dE_dx, alpha, mu, T);

    % Compute acceleration (a = F/m)
    a_y = F_de / (4/3 * pi * r_droplet^3 * rho_water) + g;  % Acceleration in y due to DEP force and gravity
    a_x = F_de / (4/3 * pi * r_droplet^3 * rho_water);     % Acceleration in x due to DEP force

    % Update velocity and position using basic Euler integration
    velocities(t, 1) = velocities(t-1, 1) + a_x * dt;
    velocities(t, 2) = velocities(t-1, 2) + a_y * dt;

    positions(t, 1) = positions(t-1, 1) + velocities(t, 1) * dt;
    positions(t, 2) = positions(t-1, 2) + velocities(t, 2) * dt;
end

% Plotting the motion of water droplet under electric field
figure(8);
plot(0:dt:t_end-dt, positions);
xlabel('Time (s)');
ylabel('Position (cm)');
title('Motion of Water Droplet Under Electric Field');
legend("x-direction","y-direction");

% Printing out the final location of the water droplet in the x and y axis
fprintf('Final position of the water droplet:\n');
fprintf('X: %f m\n', positions(end, 1));
fprintf('Y: %f m\n', positions(end, 2));

% Plotting the velocity of water droplet under electric field
figure(9);
plot(0:dt:t_end-dt, velocities);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity of Water Droplet Under Electric Field');
legend("x-direction","y-direction");

% Plotting the subplot
% Creating a new figure for the subplots
figure(10);

% Plotting the motion of water droplet under electric field
subplot(2, 1, 1); % This means 2 rows, 1 column, 1st plot
plot(0:dt:t_end-dt, positions);
xlabel('Time (s)');
ylabel('Position (cm)');
title('Motion of Water Droplet Under Electric Field');
legend("x-direction","y-direction");

% Plotting the velocity of water droplet under electric field
subplot(2, 1, 2); % This means 2 rows, 1 column, 2nd plot
plot(0:dt:t_end-dt, velocities);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity of Water Droplet Under Electric Field');
legend("x-direction","y-direction");

% figure(10);
% plot(positions(:, 2), positions(:, 1)); % X vs Y plot
% xlabel('Horizontal Position (m)');
% ylabel('Vertical Position (m)');
% title('Motion of Water Droplet Under Electric Field and Gravity');
% grid on;

% figure(11);
% plot(0:dt:t_end-dt, velocities);
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% title('Velocity of Water Droplet Under Electric Field');
% legend('Vertical Velocity', 'Horizontal Velocity');
% grid on;



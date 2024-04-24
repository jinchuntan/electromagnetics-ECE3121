% Written by Tan Jin Chun
% Last Modified: 4/9/2023

clear; close all; clc

% Constants
epsilon_0 = 8.854e-12;   % Vacuum permittivity (F/m)
rho_water = 1000;        % Density of water (kg/m^3)
r_droplet = 1e-4;        % Radius of droplet (m)
V = 120;                 % Applied voltage (V)
d = 1;                % Approximate distance from droplet to an electrode (m)

% Electric field due to one electrode
E_single = V / d;

% Droplet initial properties
position = [0, 0, 0];  % Initial position in x, y, z
velocity = [0, 0, 0];  % Initial velocity

% Time stepping loop for droplet motion
dt = 1e-3;
t_end = 1;  % Total simulation time (s)
trajectory = zeros(t_end/dt, 3);  % Store 3D positions over time

for t = 1:length(trajectory)
    % Compute net electric field at droplet position due to all 4 electrodes
    E_net = [0, 0, 0];  % Placeholder
    
    % Assuming the force is proportional to the field (this is a simplification)
    F_net = E_net * (4/3 * pi * r_droplet^3 * rho_water);
    
    % Compute acceleration (Newton's 2nd law)
    a = F_net / (4/3 * pi * r_droplet^3 * rho_water);
    
    % Update velocity and position (using basic Euler method)
    velocity = velocity + a * dt;
    position = position + velocity * dt;
    
    % Store the position
    trajectory(t, :) = position;
end

% Visualizing the droplet's motion in 3D
figure(10);

% Plot trajectory
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'b');
hold on;

% Plotting the water stream as a blue line from the top to bottom
z_start = 0;  % Top of the Z-axis (above the electrode)
z_end = -5;   % Bottom of the Z-axis
x_water = 0;  % X-coordinate for the water stream, assuming it's at the center
y_water = 0;  % Y-coordinate for the water stream, assuming it's at the center
plot3([x_water, x_water], [y_water, y_water], [z_start, z_end], 'b-', 'LineWidth', 2);  % Water stream

% Plot electrodes
% Initialising the water droplet movement
plot3([-d, d, 0, 0], [0, 0, -d, d], [0, 0, 0, 0], 'rs', 'MarkerSize', 5, 'LineWidth', 2);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Trajectory of Water Droplet with Water Stream');
grid on;
axis equal;
hold off

% View from top
figure(11);
plot(trajectory(:, 1), trajectory(:, 2), 'b');
hold on;
plot([-d, d, 0, 0], [0, 0, -d, d], 'rs', 'MarkerSize', 5, 'LineWidth', 2);

% Adding a blue circle at the center of the four electrodes
circle_center = [0, 0];
circle_radius = (d/10)/2;  % Half the distance d
viscircles(circle_center, circle_radius, 'LineWidth', 1, 'EdgeColor', 'b');

% Labelling the plot
xlabel('X');
ylabel('Y');
title('Top View of Droplet Trajectory');
axis equal;



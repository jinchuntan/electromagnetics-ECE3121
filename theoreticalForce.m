function F_elec = theoreticalForce(E, dE_dx, alpha, mu, T)
    % Constants
    k = 1.38e-23; % Boltzmann constant in J/K

    % Calculating the electric force based on the given equation
    F_elec = - (alpha + (mu^2 / (3 * k * T))) * E * dE_dx;
end
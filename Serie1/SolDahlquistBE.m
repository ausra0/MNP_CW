function [solutions, temps] = SolDahlquistBE(lamb, y0, T, N)
% lamb is a parameter of the ODE (can be complex)
% y0 is the initial point 
% [0, T] interval of time on which we are solving the ODE 
% N the number of time steps 

    % definir les paramètres pour la résolution 
    h = T/N; 
    temps = linspace(0, T, N+1); 
    
    % définir un vecteur pour les solutions
    solutions = zeros(1, N+1); 
    
    % se souvenir de la condition initiale 
    y = y0; 
    solutions(1) = y0; 
    
    % calculer la solution 
    for i = 1:N
       y = 1/(1 - lamb*h)*y; 
       solutions(i+1) = y; 
    end
end
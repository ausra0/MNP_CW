function [tu, u, u_analyt] = solveHeat(L, J, T, N, nu, u, method)
    % résoud l'équation de la chaleur numériquement et analytiquement sur
    % [0, L]x[0, T]
    % -----------------------------------------------------------
    % [0, L] intervalle d'espace considéré
    % J nombre de points considérés 
    % [0, T] intervalle de temps considéré
    % N nombre de points considérés
    % nu constante 
    % u consition initiale (function handle)
    % method : 'FE' = forward euler/ 'BE' = backward euler 
    % -----------------------------------------------------------
    
    % définir les paramètres internes 
    Nstep = J;
    
    % construct A 
    A = construireA(L, J, nu); 
    
    % u0 condition initiale u(0, t)
    u0 = arrayfun(u, linspace(0, L, J)); 
    u0 = u0'; 
    
    
    % résoudre
    switch method 
        case 'FE'
            [tu, u] = ForwardEulerLin(A, u0, T, N);
        case 'BE'
            [tu, u] = BackwardEulerLin(A, u0, T, N);
    end
    
    % calculer la solution numérique 
    x = linspace(0, L, J); 
    t = linspace(0, T, N); 
    u_analyt = calcSolAnalytique(Nstep, x, t, nu, L);
end
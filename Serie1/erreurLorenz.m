function erreurLorenz(sig, rho, beta, x0, y0, z0, delt, T, Nexact)
    % calcule les erreurs de troncature locale et globale. 
    % sig, rho, beta sont les paramètres du système de Lorenz 
    % [x0, y0, z0]' la condition initiale du système 
    % delt les incréments en temps pour lesquels on calcule l'erreur 
    % [0, T] l'échelle de temps considéré
    % Nexact pas de temps pour obtenir une solution approchée à la solution
    %   exacte
        
    % définir un tableau pour contenir les erreures ... 
    % ... locales 
    err_loc = zeros(1, length(delt)); 
    % ... globales 
    err_glob = zeros(1, length(delt)); 
    
    % calculer la solution exacte 
    sol_exact = resoudreLorenz(sig, rho, beta, x0, y0, z0, T, Nexact, 0); 
    
    % calculer les erreurs locales et globales
    for dt=1:length(delt)
        % pas de temps pour chacune des méthodes
        N = T/delt(dt); 
        
        % calculer la solution approchée
        sol = resoudreLorenz(sig, rho, beta, x0, y0, z0, T, N, 0); 
        
        % calculer l'erreur globale 
        err_glob(dt) = norm(sol(:, end) - sol_exact(:, end), 2); 
        
        % calculer la solution exacte pour un pas de la méthode considérée 
        t1 = delt(dt); % un pas 
        Nt1 = floor(Nexact/T); % nombre de pas 
        sol_ex = resoudreLorenz(sig, rho, beta, x0, y0, z0, t1, Nt1, 0); 
        
        % calculer l'erreur locale 
        err_loc(dt) = norm(sol(:, 2) - sol_ex(:, end), 2); 
    end
    
    figure
    % afficher l'erreur locale 
    subplot(1, 2, 1); 
    loglog(delt, err_loc, '-ob'); 
    hold on; 
    loglog(delt, delt.^2, '--r')
    hold off; 
    legend('Erreur locale', 'O(h^{2})'); 
    
    % afficher l'erreur globale
    subplot(1, 2, 2); 
    loglog(delt, err_glob, '-ob'); 
    hold on; 
    loglog(delt, delt, '--r')
    hold off; 
    legend('Erreur globale', 'O(h)'); 
end
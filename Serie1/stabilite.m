function stabilite(sig, rho, beta, x, y, z)
    % calcule si le point (x, y, z) est un point fixe stable pour le
    % système de Lorenz 
    % sig, rho, beta sont les paramètres de l'EDO 
    
    Jac = [-sig, sig, 0; rho, -1, -x; y, x, -beta]; 

    E = eig(Jac); 

    if(prod(E<0)==1)
        fprintf('Le système est stable pour le point (%.2f, %.2f, %.2f)\n', x, y, z); 
    else 
        fprintf('Le système est instable pour le point (%.2f, %.2f, %.2f)\n', x, y, z); % really ? 
    end
end
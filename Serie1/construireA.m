function A = construireA(L, J, nu)
    % construit la matrice A
    % [0, L] le segment d'espace sur lequel on défini notre équation 
    % J le nombre de points sur la discretisation
    
    % pas 
    h = L/(J+1);
    
    % construire la matrice
    e = ones(J, 1); 
    A = spdiags([e, -2*e, e], -1:1, J, J);
    
    % ajouter les coefficients
    A = nu/(h^2).*A; 
end
function u = calcSolAnalytique(Nf, x, t, nu, L)
% calcule la solution analytique pour l'équation de la chaleur et u0 = 20
% (tronquée aux k premiers termes)
% x : un vecteur de valeurs d'espace prises (vect. ligne)
% t : un vecteur de valeurs de temps prises (vect. ligne)

% calculer le nombre pas à faire
n = ceil(Nf/2);

% initialiser la solution:
u = zeros(length(x), length(t)); 
for k=1:n 
    sig = pi*(2*k-1)/L; % 2k - 1 = j car j est impair
    a = (80/(sig*L)); 
    b = exp((-sig^2/nu).*t); 
    c = sin(sig.*x'); 
    u = u + a.* c.*b; 
end

end
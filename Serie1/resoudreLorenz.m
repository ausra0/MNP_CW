function solution = resoudreLorenz(sig, rho, beta, x0, y0, z0, T, N, disp_switch)
    % Resoud le système de Lorenz à l'aide de la méthode d'euler explicite
    % sig, rho, beta sont les paramètres de l'EDO (scalaires) 
    % (x0, y0 , z0), le point initial
    % [0, T] l'intervalle de temps sur lequel on résoud le système 
    % N est le pas de temps
    % disp_switch boolean for displaying or not the figure
    
    % définir les paramètres pour la résolution 
    h = T/N; 
    x = x0; 
    y = y0; 
    z = z0; 
    
    % garder la solution en mémoire
    solution = zeros(3, N+1); 
    solution(:, 1) = [x0; y0; z0]; 
    
    % résoudre de manière itérative 
    for i = 1:N
        % résoudre
        x = x + h*sig*(y - x); 
        y = y + h*(x*(rho - z)- y); 
        z = z + h*x*y - h*beta*z; 
        
        % garder la solution en mémoire
        solution(:,i+1) = [x; y; z]; 
    end 
    
    % afficher solution 
    if(disp_switch==1)
        figure(1)
        plot3(solution(1, :), solution(2, :), solution(3, :)); 
        title('Trajectoire pour x(t) (système de Lorenz)');
    end
end
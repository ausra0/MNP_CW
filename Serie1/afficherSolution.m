function afficherSolution(solution, N, T)
    % solution est une matrice 3xN
    
    x = linspace(0, T, N+1); 
    
    figure 
    title('Différence des trajectoires (système de Lorenz)');
    
    subplot(3, 1, 1); 
    plot(x, solution(1, :)); 
    xlabel('temps'); 
    ylabel('|x_{1}(t) - x_{2}(t)|'); 
    
    subplot(3, 1, 2); 
    plot(x, solution(2, :));
    xlabel('temps'); 
    ylabel('|y_{1}(t) - y_{2}(t)|'); 
    
    subplot(3, 1, 3); 
    plot(x, solution(3, :));
    xlabel('temps'); 
    ylabel('|z_{1}(t) - z_{2}(t)|'); 
    
end
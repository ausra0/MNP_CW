function afficherStabFE(dt, lambi)
% dt is the time step 

    figure
    drawCircle(-1/dt, 0, 1/dt); 
    
    
    hold on; 
    % draw the x-y axis
    plot([-5/dt, 5/dt], [0, 0], 'k-');
    plot([0, 0], [-5/dt, 5/dt], 'k-'); 
    
    % draw the lines
    col = ['m', 'g', 'c', 'b', 'y'];
    for j=1:length(lambi)
        plot([0, 5*real(lambi(j))/(dt*abs(lambi(j)))], [0, 5*imag(lambi(j))/(dt*abs(lambi(j)))], col(j))
    end
    hold off;
    
    title('Zone de stabilité en rouge pour FE');
    axis equal
    grid on
    
end

function drawCircle(x, y, r)
    % draw the circle
    rectangle('Position',[x-r, y-r, 2*r, 2*r] ,'Curvature',[1 1], 'FaceColor', 'red', 'EdgeColor', 'red'); 
     
end
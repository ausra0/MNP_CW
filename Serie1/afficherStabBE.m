function afficherStabBE(dt)
% dt is the time step 

    figure
    set(gca,'Color','red')
    drawCircle(1/dt, 0, 1/dt); 

    % draw the x-y axis
    hold on; 
    plot([-5/dt, 5/dt], [0, 0], 'k-');
    plot([0, 0], [-5/dt, 5/dt], 'k-'); 
    hold off; 
    
    title('Zone de stabilité en rouge pour BE');
    axis equal
    grid on
    
end

function drawCircle(x, y, r)
    % draw the circle
    rectangle('Position',[x-r, y-r, 2*r, 2*r] ,'Curvature',[1 1], 'FaceColor', 'white', 'EdgeColor', 'white'); 
   
end
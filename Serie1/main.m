% Serie 1 - Méthodes Numériques pour le Parallélisme en Temps 
% Ausra Pogozelskyte 
%% SYSTEME DE LORENZ 

%% Stabilité des points fixes

% Définir le système 
sig = 10; 
rho = 28; 
beta = 8/3; 

% premier point fixe 
x = sqrt(beta*(rho - 1)); 
y = x; 
z = rho - 1; 

stabilite(sig, rho, beta, x, y, z); 

% deuxième point fixe 
x = -sqrt(beta*(rho - 1)); 
y = x; 
z = rho - 1; 

stabilite(sig, rho, beta, x, y, z); 

% on remarque que les deus points fixes ont la même stabilité comme
% démontré au point (d)

%% Resoudre Lorenz via Euler Implicite 

% Définir le système 
sig = 10; 
rho = 28; 
beta = 8/3; 

% Définit les points fixes 
xf1 = 0; 
yf1 = 0; 
zf1 = 0; 

xf2 = sqrt(beta*(rho - 1)); 
yf2 = x; 
zf2 = rho - 1; 

xf3 = -sqrt(beta*(rho - 1)); 
yf3 = x; 
zf3 = rho - 1;

% point initial 
x0 = 20; 
y0 = 5; 
z0 = -5; 

% paramètres pour la résolution 
T= 20; 
N = 20000; 
epsi = 1e-3; % modification des conditions initiales % 6e-10

solution1 = resoudreLorenz(sig, rho, beta, x0, y0, z0, T, N, 1); 

hold on; 
solution2 = resoudreLorenz(sig, rho, beta, x0 + epsi, y0, z0, T, N, 1);
plot3([xf1, xf2, xf3], [yf1, yf2, yf3], [zf1, zf2, zf3], '.k'); 
legend('solution', 'solution modifiée', 'points fixes'); 
hold off; 

%% Afficher les différences entre les solutions 
sol = abs(solution1 - solution2); 
afficherSolution(sol, N, T); 

%% Calcul de l'erreur de troncateur locale 
% Définir le système 
sig = 10; 
rho = 28; 
beta = 8/3;

% point initial 
x0 = 20; 
y0 = 5; 
z0 = -5; 

% Initialiser les paramètres de résolution
delt = [1e-4, 1e-3, 1e-2, 1e-1]; %dt = T/Nstep les incréments de temps pour lequels on va tester 
T = 0.1; 
Nexact = 1e5; 

% calculer les erreurs locales et globales
erreurLorenz(sig, rho, beta, x0, y0, z0, delt, T, Nexact)

%% EQUATION DE DAHLQUIST 

%% Région de stabilité Forward Euler
dt = T/N; %1; 
lambi = [-1, 1i, 1i - 2]; 
afficherStabFE(dt, lambi); 

%% Région de stabilité Backward Euler
dt = T/N; 
afficherStabBE(dt); 

%% Solution de l'équation de Dahlquist 
 
% définir les paramètres de notre problème 
lamb = [1i, 1i - 1]; 
y0 = 1; 
T = 1; 
N = 10; 

for j= 1:length(lamb)
    % Resoudre eq. Dahlquist via Forward Euler 
    [sol1, temps] = SolDahlquistFE(lamb(j), y0, T, N); 

    % Résoudre eq. Dahlquist via Backward Euler
    [sol2, ~] = SolDahlquistBE(lamb(j), y0, T, N); 
    
    % calculer la solution réelle
    f = @(t)(exp(lamb(j).*t));
    solre = f(temps); 
    
    
    % afficher les solutions dans le plan complexe
    figure 
    plot(real(sol1), imag(sol1));
    hold on
    plot(real(sol2), imag(sol2)); 
    plot(real(solre), imag(solre));
    hold off
    
    % mise en forme pour le plot précédent
    title(sprintf('Solutions de l''équation de Dahlquist pour lambda = %d + %d i', real(lamb(j)), imag(lamb(j)))); 
    legend('Forward Euler', 'Backward Euler', 'real solution'); 
    xlabel('Re(y)'); 
    ylabel('Im(y)'); 
    
    
    % afficher les modules et phases des solutions 
    
    % afficher d'abord le module de la solution en fonction du temps
    figure 
    subplot(2, 1, 1); 
    plot(temps, abs(sol1)); 
    hold on 
    plot(temps, abs(sol2)); 
    plot(temps, abs(solre)); 
    hold off 
    
    % mise en forme pour le plot précédent
    title(sprintf('Modules de la sol. de l''équation de Dahlquist pour lambda = %d + %d i', real(lamb(j)), imag(lamb(j)))); 
    legend('Forward Euler', 'Backward Euler', 'real solution'); 
    xlabel('t'); 
    ylabel('|y(t)|'); 
    
    % afficher la phase de la solution en fonction du temps 
    subplot(2, 1, 2); 
    plot(temps, angle(sol1));
    hold on 
    plot(temps, angle(sol2)); 
    plot(temps, angle(solre)); 
    hold off
    
    % mise en forme pour le plot précédent 
    title(sprintf('Phases de la sol. de l''équation de Dahlquist pour lambda = %d + %d i', real(lamb(j)), imag(lamb(j)))); 
    legend('Forward Euler', 'Backward Euler', 'real solution'); 
    xlabel('t'); 
    ylabel('\theta(y(t)) in radians'); 
    
    
    % afficher les erreurs de modules et de phases des solutions 
    
    % afficher d'abord l'erreur de module de la solution en fonction du temps
    figure 
    subplot(2, 1, 1); 
    plot(temps, abs(abs(sol1) - abs(solre))); 
    hold on 
    plot(temps, abs(abs(sol2) - abs(solre))); 
    hold off 
    
    % mise en forme pour le plot précédent
    title(sprintf('Erreurs de modules pour lambda = %d + %d i', real(lamb(j)), imag(lamb(j)))); 
    legend('Forward Euler', 'Backward Euler'); 
    xlabel('t'); 
    ylabel('||y(t)| - |y_{approx}(t)||'); 
    
    % afficher la phase de la solution en fonction du temps 
    subplot(2, 1, 2); 
    plot(temps, abs(angle(sol1) - angle(solre)));
    hold on 
    plot(temps, abs(angle(sol2) - angle(solre)));  
    hold off
    
    % mise en forme pour le plot précédent 
    title(sprintf('Erreur de phases pour lambda = %d + %d i', real(lamb(j)), imag(lamb(j)))); 
    legend('Forward Euler', 'Backward Euler'); 
    xlabel('t'); 
    ylabel('|\theta(y(t)) -\theta(y_{approx}(t))| in radians'); 
    
end

%% EQUATION DE LA CHALEUR 

%% Solution exacte
% définir les paramètres du problème
L = 1; 
T = 1/2; 
x = linspace(0, L, 50); 
t = linspace(0, T, 20); 
nu = 1; 

% définir les rangs pour lesquels on veut calculer la solution analytique
Nf = [5, 10, 20, 40]; % 100

for k = Nf
    % calculer la solution 
    u = calcSolAnalytique(k, x, t, nu, L); 
    
    % afficher la solution 
    figure 
    contour(t, x, u); 
    title(sprintf('Solution analytique jusqu''au rang Nf = %d', k)); 
    ylabel('space'); 
    xlabel('time'); 
    
end

%% La matrice A et ses valeurs propres

% paramètres du problème 
L = 11; 
J = 10;
nu = 1; % mu est en général né

% construire la matrice A 
A = construireA(L, J, nu);

% calculer les valeurs propres de A 
E = eig(A); 
disp(E); 

% calculer les valeurs propres théoriques de A 
h = L/(J+1); 
cst = nu/h; 
a = -2*cst; 
b = 1*cst; 
c = b; 
Et = a + 2*sqrt(b*c).*cos((1:J).*pi./(J+1));

% calculer la différence entre les valeurs propres théoriques et empiriques
disp(abs(sort(Et)' - E)>1e-14); 

%% Forward Euler : Résolution numérique de l'éq. de la chaleur

% Paramètres de notre problème 
L = 1; 
J = 99;
T = 1/2; 
N = 1e4;
u0 = @(x)(20); 
nu = 1; 
method = 'FE'; 

% Afficher dt
fprintf('df = %0.5f\n', T/N); 

% résoudre l'équation de la chaleur 
[t, u, u_analyt] = solveHeat(L, J, T, N, nu, u0, method); 

% plot the output
figure 
subplot(1, 2, 1); 
contour(linspace(0, T, N+1), linspace(0, L, J), u); 
xlabel('time'); 
ylabel('space'); 
title('Numerical solution of Heat Equation (FE)'); 

subplot(1, 2, 2); 
contour(linspace(0, T, N), linspace(0, L, J), u_analyt); 
xlabel('time'); 
ylabel('space'); 
title('Analytic solution of Heat Equation (FE)'); 

% afficher la différence entre les solutions 
figure 
prec = 8; 
contour(linspace(0, T, N), linspace(0, L, J), abs(u(:, 2:end) - u_analyt), logspace(-prec, 0, prec +1));
xlabel('time'); 
ylabel('space'); 
title('Error between solutions')

%% Forward Euler : Nstep réduit

% Paramètres de notre problème 
L = 1; 
J = 99;
T = 1/2; 
N = 1e4 - 6;
u0 = @(x)(20); 
nu = 1; 
method = 'FE'; 

% Afficher dt
fprintf('df = %0.5f\n', T/N); 

% résoudre l'équation de la chaleur 
[~, u, u_analyt] = solveHeat(L, J, T, N, nu, u0, method); 

% plot the output
figure 
subplot(1, 2, 1); 
contour(linspace(0, T, N+1), linspace(0, L, J), u); 
xlabel('time'); 
ylabel('space'); 
title('Numerical solution of Heat Equation (FE)'); 

subplot(1, 2, 2); 
contour(linspace(0, T, N), linspace(0, L, J), u_analyt); 
xlabel('time'); 
ylabel('space'); 
title('Analytic solution of Heat Equation (FE)'); 

% afficher la différence entre les solutions 
figure 
prec = 8; 
contour(linspace(0, T, N), linspace(0, L, J), abs(u(:, 2:end) - u_analyt), logspace(-prec, 0, prec +1));
xlabel('time'); 
ylabel('space'); 
title('Error between solutions')

%% Forward Euler : J réduit

% Paramètres de notre problème 
L = 1; 
J = 49;
T = 1/2; 
N = 2.5e3;
u0 = @(x)(20); 
nu = 1; 
method = 'FE'; 

% Afficher dt
fprintf('df = %0.5f\n', T/N); 

% résoudre l'équation de la chaleur 
[~, u, u_analyt] = solveHeat(L, J, T, N, nu, u0, method); 

% plot the output
figure 
subplot(1, 2, 1); 
contour(linspace(0, T, N+1), linspace(0, L, J), u); 
xlabel('time'); 
ylabel('space'); 
title('Numerical solution of Heat Equation (FE)'); 

subplot(1, 2, 2); 
contour(linspace(0, T, N), linspace(0, L, J), u_analyt); 
xlabel('time'); 
ylabel('space'); 
title('Analytic solution of Heat Equation (FE)'); 

% afficher la différence entre les solutions 
figure 
prec = 8; 
contour(linspace(0, T, N), linspace(0, L, J), abs(u(:, 2:end) - u_analyt), logspace(-prec, 0, prec +1));
xlabel('time'); 
ylabel('space'); 
title('Error between solutions')

%% Backward Euler : Résolution numérique de l'éq. de la chaleur

% Paramètres de notre problème 
L = 1; 
J = 99;
T = 1/2; 
N = 4e0; % 5e3
u0 = @(x)(20); 
nu = 1; 
method = 'BE'; 

% résoudre l'équation de la chaleur 
[~, u, u_analyt] = solveHeat(L, J, T, N, nu, u0, method); 

% plot the output
figure 
subplot(1, 2, 1); 
contour(linspace(0, T, N+1), linspace(0, L, J), u); 
xlabel('time'); 
ylabel('space'); 
title('Numerical solution of Heat Equation (BE)'); 

subplot(1, 2, 2); 
contour(linspace(0, T, N), linspace(0, L, J), u_analyt); 
xlabel('time'); 
ylabel('space'); 
title('Analytic solution of Heat Equation (BE)'); 

% afficher la différence entre les solutions 
figure 
prec = 8; 
contour(linspace(0, T, N), linspace(0, L, J), abs(u(:, 2:end) - u_analyt), logspace(-prec, 0, prec +1));
xlabel('time'); 
ylabel('space'); 
title('Error between solutions')
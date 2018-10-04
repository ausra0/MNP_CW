function [t, u] = ForwardEulerLin(A, u0, T, N)
%FORWARDEULERLIN Solve a linear system of ODE using Forward Euler
%    Solve a linear system of ODE using Forward Euler.
%    Considering the problem
%
%    .. math::
%        \\frac{dU}{dt} = AU,
%
%    computes the numerical solution with Forward Euler, between :math:`t=0`
%    and **T**, with **u0** the initial solution, using **N** steps.
%
%    Parameters
%    ----------
%    A : matrix of size JxJ
%        The matrix of the linear system
%    u0 : vector of size J
%        The initial vector
%    T : float
%        The final time of the solution
%    N : int
%        The number of numerical time steps
% 
%     Returns
%     -------
%     u : matrix of size JxN
%         The solution at each time steps (including initial solution)
%     t : vector of size N
%         The times of the solutions

% pas de temps
dt = T/N;
J = length(u0);
% solution numérique (en fnct du temps)
u = zeros(J, N+1);
% condition initiale pour le temps
u(:, 1) = u0;
% prend en compte la partie temporelle
R = speye(J)+dt*A;
for i=1:N
    u(:, i+1) = R*u(:, i);
end
t = linspace(0, T, N+1);
end

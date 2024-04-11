clc; 
clear; 
close all;

import Skeleton_NS_solver.*

% solve_flow(N, Re, tol, K=1, CFL) % standard Re = 1000, tol = 10^-4
% Outputs: u, p, xi, diff, iter, error

%% Timestep/Tolerance analysis
% CFLs = [0.1, 0.3, 0.5, 0.7, 0.9];
% Ks = [1, 2, 3, 4, 5];
% tols = [10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8];
% 
% elapsed_times = zeros(1, length(tols));
% errors = zeros(1, length(tols));
% 
% for i = 1:length(tols)
%     tol = tols(i)
%     tic;
%     results = Skeleton_NS_solver(45, 1000, tol, 1);
%     elapsed_time = toc;
%     elapsed_times(i) = elapsed_time;
%     errors(i) = results.error;
% end

% % Plot elapsed_times against CFLs
% hold on
% % plot(tols, elapsed_time0*ones(1, length(Ks)));
% semilogx(tols, elapsed_times, 'o-');
% xlabel('Tolerance');
% ylabel('Elapsed Time (seconds)');
% title('Elapsed Time vs. Tolerance');
% hold off
% 
% hold on
% % plot(Ks, error0*ones(1, length(Ks)));
% semilogx(tols, errors, 'o-');
% xlabel('Tolerance');
% ylabel('Errors ');
% title('Errors(RMS) vs. Tolerance');
% hold off

%% tolerance visual verification
% N = 45;
% result = Skeleton_NS_solver(N, 1000, 10^-5, 1);
% [Ht11, H1t1] = hodges11(result.h, result.th, N);
% [tE10] = maketE10(N);
% u_prim = Ht11*result.u;
% psi = linsolve(full(tE10),u_prim);
% psistag = rot90(flipud(reshape(psi, N+1, N+1)), -1); % Streamfunction
% 
% % Stagger mesh and Reshape for plotting
% [pX,pY] = staggered(result.x); % Primal
% 
% levels3 = [0.1175 0.115 0.11 0.1 0.09 0.07 0.05 0.03 0.01 10^(-4) 10^(-5) 10^(-10) 0 ...
%     -10^(-6) -10^(-5) -5*10^(-5) -10^(-4) -2.5*10^(-4) -5*10^(-4) -10^(-3) -1.5*10^(-3)];
% title3 = ['Streamfunction Field (N = 45, tol = 10^-5)'];
% contour(pX, pY, psistag, 'LevelList', levels3);
% title(title3);


%% Default dt vs Courant dt
% [u, p, xi, diff, iter, error] = Skeleton_NS_solver(47, 1000, 10^-5, 0.7);
% error_courant = error;
% [u, p, xi, diff, iter, error] = Skeleton_NS_solver(47, 1000, 10^-5);
% error_default = error;
% % Resize lists to have the same length (max length_
% max_length = max(length(diff_listCFL), length(diff_list));
% diff_listCFL = [diff_listCFL, zeros(1, max_length - length(diff_listCFL))];
% diff_list = [diff_list, zeros(1, max_length - length(diff_list))];
% % Create x-axis vector
% x = 1:max_length;


%% Pressure matrix A
% result = Skeleton_NS_solver(55, 1000, 10^-5, 1);
% 
% A = result.A; % Pressure matrix
% 
% disp(['Row sum is zero: ', num2str(all(sum(A, 2) == 0))]); % Check if row sum is zero
% disp(['Matrix is symmetric: ', num2str(isequal(A, A'))]); % Check if symmetric
% disp(['Matrix is singular: ', num2str((rank(full(A)) < min(size(A))))]); % Check if singular
% 
% [V, D] = eigs(A'*A, 1);
% eigenmode = A * V;
% eigvector = eigenmode / norm(eigenmode);


%% Integrated vorticity
N = 35;
Res = [200, 600, 1000, 1400, 1800];
tol = 10^-3;
intvort = zeros(1, length(Res));

for i = 1:length(Res)
    Re = Res(i)
    results = Skeleton_NS_solver(N, Re, tol, 1);
    xi_other = results.xiother;
    intvort(i) = sum(abs(xi_other)); % Is equal to 1
end

plot(Res, intvort, 'o-');
xlabel('Re [-]');
ylabel('Integrated vorticity [m*m/s]');
title('Integrated vorticity vs. Reynolnds number');


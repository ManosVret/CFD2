clc; 
clear; 
close all;
filename='results_N_63.mat'; % Ν = 15, 31, 47, 55, 63
load(filename)

%% Preprocess
% Remove boundary points of pressure
p([(N+2)^2-4-N+1:(N+2)^2-4]) = [];
p([1:N]) = [];
for i=1:N
    p([1+(i-1)*(N+2)]) = 9;
    p([1+(i-1)*(N+2)+(N+1)]) = 9;
end
p = p(p ~= 9);

% Create mesh with segment widths (dual)
hmesh = zeros(size(u));
for i=1:N
    hmesh([(i-1)*(N+1)+1:(i-1)*(N+1)+1+N]) = h;
end
for i=1:N+1
    hmesh([N*(N+1)+1+(i-1)*N:N*(N+1)+1+(i-1)*N+N-1]) = h(i);
end

% Velocity interpolation on the pressure points (dual grid)
umesh = u./hmesh;
uinterp = zeros(size(p));
uxinterp = zeros(size(p));
uyinterp = zeros(size(p));
for i=1:N
    for j=1:N  
        ux = umesh((i-1)*(N+1) + j ) + (umesh((i-1)*(N+1) + j + 1)-umesh((i-1)*(N+1) + j)) ... 
            /(hmesh((i-1)*(N+1) + j)+hmesh((i-1)*(N+1) + j + 1)) * hmesh((i-1)*(N+1) + j);
        uxinterp((i-1)*N + j) = ux;
        uy = umesh((i-1)*N + j + N*(N+1)) + (umesh((i-1)*N + j + N*(N+1) + N)-umesh((i-1)*N + j + N*(N+1))) ...
            /(hmesh((i-1)*N + j + N*(N+1))+hmesh((i-1)*N + j + N*(N+1) + N)) * hmesh((i-1)*N + j + N*(N+1));
        uyinterp((i-1)*N + j) = uy;
        uinterp((i-1)*N + j) = sqrt(ux.^2+uy.^2);
    end
end

% Remove dynamic pressure and normalize for center point
pdyn = p - 0.5*uinterp.^2;
pdyn = pdyn - pdyn(round((N^2+1)/2));

% Streamfunction (tE10*psi = u_prim -> solve for psi)
u_prim = Ht11*u;
psi = linsolve(full(tE10),u_prim);

% Stagger mesh and Reshape for plotting
[dX,dY] = staggered(tx); % Dual
[pX,pY] = staggered(x); % Primal

ustag = rot90(flipud(reshape(uinterp,N,N)), -1); % Velocity
uxstag = rot90(flipud(reshape(uxinterp,N,N)), -1); % Velocity x
uystag = rot90(flipud(reshape(uyinterp,N,N)), -1); % Velocity y
pstag = rot90(flipud(reshape(pdyn, N, N)), -1); % Pressure
xistag = rot90(flipud(reshape(xi,N+1,N+1)), -1); % Vorticity
psistag = rot90(flipud(reshape(psi, N+1, N+1)), -1); % Streamfunction

%% Plotting the Fields
% Define the contour levels and titles for each plot
levels1 = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0 -0.002];
title1 = ['Pressure Field (N = ', num2str(N), ')'];

levels2 = [5 4 3 2 1 0.5 0 -0.5 -1 -2 -3];
title2 = ['Vorticity Field (N = ', num2str(N), ')'];

levels3 = [0.1175 0.115 0.11 0.1 0.09 0.07 0.05 0.03 0.01 10^(-4) 10^(-5) 10^(-10) 0 ...
    -10^(-6) -10^(-5) -5*10^(-5) -10^(-4) -2.5*10^(-4) -5*10^(-4) -10^(-3) -1.5*10^(-3)];
title3 = ['Streamfunction Field (N = ', num2str(N), ')'];

% Create a new figure for each contour plot and save it in the "Plots" folder
plotsFolderPath = fullfile(fileparts(mfilename('fullpath')), 'Plots');

% Save the plots in the "Plots" folder
filename1 = fullfile(plotsFolderPath, ['pressure_field_N', num2str(N), '.png']);
figure(1);
contour(dX, dY, pstag, 'LevelList', levels1);
title(title1);
saveas(gcf, filename1);

filename2 = fullfile(plotsFolderPath, ['vorticity_field_N', num2str(N), '.png']);
figure(2);
contour(pX, pY, xistag, 'LevelList', levels2);
title(title2);
saveas(gcf, filename2);

filename3 = fullfile(plotsFolderPath, ['streamfunction_field_N', num2str(N), '.png']);
figure(3);
contour(pX, pY, psistag, 'LevelList', levels3);
title(title3);
saveas(gcf, filename3);

%% Plotting the Halflines
yref = flip([1.00000 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5000 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0.0000]);
uref = flip([-1.0000000 -0.6644227 -0.5808359 -0.5169277 -0.4723329 -0.3372212 -0.1886747 -0.0570178 0.0620561 0.1081999 0.2803696 0.3885691 0.3004561 0.2228955 0.2023300 0.1812881 0.0000000]);
p_vary_y = flip([0.052987 0.052009 0.051514 0.050949 0.050329 0.034910 0.012122 -0.000827 0.000000 0.004434 0.040377 0.081925 0.104187 0.108566 0.109200 0.109689 0.110591]);
vort_vary_y = flip([14.7534 12.0670 9.49496 6.95968 4.85754 1.76200 2.09121 2.06539 2.06722 2.06215 2.26772 1.05467 -1.63436 -2.20175 -2.31786 -2.44960 -4.16648]);

xref = [0.0000 0.0312 0.0391 0.0469 0.0547 0.0937 0.1406 0.1953 0.5000 0.7656 0.7734 0.8437 0.9062 0.9219 0.9297 0.9375 1.0000];
vref = [0.0000000 -0.2279225 -0.2936869 -0.3553213 -0.4103754 -0.5264392 -0.4264545 -0.3202137 0.0257995 0.3253592 0.3339924 0.3769189 0.3330442 0.3099097 0.2962703 0.2807056 0.0000000];
p_vary_x = [0.077455 0.078837 0.078685 0.078148 0.077154 0.065816 0.049029 0.034552 0.000000 0.044848 0.047260 0.069511 0.084386 0.086716 0.087653 0.088445 0.090477];
vort_vary_x = [-5.46217 -8.44350 -8.24616 -7.58524 -6.50867 0.92291 3.43016 2.21171 2.06722 2.06122 2.00174 0.74207 -0.82398 -1.23991 -1.50306 -1.83308 -7.66369];

% Plots at x = 0.5
filename4 = fullfile(plotsFolderPath, ['pressure_Xmid_N', num2str(N), '.png']);
figure(4);
plot(cumsum(th), pstag(:, floor(N/2))); % Pressure
hold on; 
plot(yref, p_vary_y, 'rs', 'MarkerFaceColor', 'red'); % Reference values
hold off; 
title(['Pressure at X=0.5 (N = ', num2str(N), ')']);
legend('Calculated', 'Reference Values', 'Location', 'best');
saveas(gcf, filename4);

filename5 = fullfile(plotsFolderPath, ['vorticity_Xmid_N', num2str(N), '.png']);
figure(5);
plot(cumsum(h), xistag(:, floor(N/2))); % Vorticity
hold on; 
plot(yref, vort_vary_y, 'rs', 'MarkerFaceColor', 'red'); % Reference values
hold off; 
title(['Vorticity at X=0.5 (N = ', num2str(N), ')']);
legend('Calculated', 'Reference Values', 'Location', 'best');
saveas(gcf, filename5);

filename6 = fullfile(plotsFolderPath, ['ux_Xmid_N', num2str(N), '.png']);
figure(6);
plot(cumsum(th), uxstag(:, floor(N/2))); % Velocity x
hold on; 
plot(yref, uref, 'rs', 'MarkerFaceColor', 'red'); % Reference values
hold off; 
title(['Horizontal velocity at X=0.5 (N = ', num2str(N), ')']);
legend('Calculated', 'Reference Values', 'Location', 'best');
saveas(gcf, filename6);

filename7 = fullfile(plotsFolderPath, ['uy_Xmid_N', num2str(N), '.png']);
figure(7);
plot(cumsum(th), uystag(:, floor(N/2))); % Velocity y
title(['Vertical velocity at X=0.5 (N = ', num2str(N), ')']);
saveas(gcf, filename7);

% Plots at y = 0.5
filename8 = fullfile(plotsFolderPath, ['pressure_Ymid_N', num2str(N), '.png']);
figure(8);
plot(cumsum(th), pstag(floor(N/2), :)); % Pressure
hold on; 
plot(xref, p_vary_x, 'rs', 'MarkerFaceColor', 'red'); % Reference values
hold off; 
title(['Pressure at Y=0.5 (N = ', num2str(N), ')']);
legend('Calculated', 'Reference Values', 'Location', 'best');
saveas(gcf, filename8);

filename9 = fullfile(plotsFolderPath, ['vorticity_Ymid_N', num2str(N), '.png']);
figure(9);
plot(cumsum(h), xistag(floor(N/2), :)); % Vorticity
hold on; 
plot(xref, vort_vary_x, 'rs', 'MarkerFaceColor', 'red'); % Reference values
hold off; 
title(['Vorticity at Y=0.5 (N = ', num2str(N), ')']);
legend('Calculated', 'Reference Values', 'Location', 'best');
saveas(gcf, filename9);

filename10 = fullfile(plotsFolderPath, ['ux_Ymid_N', num2str(N), '.png']);
figure(10);
plot(cumsum(th), uxstag(floor(N/2), :)); % Velocity x
title(['Horizontal velocity at Y=0.5 (N = ', num2str(N), ')']);
saveas(gcf, filename10);

filename11 = fullfile(plotsFolderPath, ['uy_Ymid_N', num2str(N), '.png']);
figure(11);
plot(cumsum(th), uystag(floor(N/2),:)); % Velocity y
hold on; 
plot(xref, vref, 'rs', 'MarkerFaceColor', 'red'); % Reference values
hold off; 
title(['Vertical velocity at Y=0.5 (N = ', num2str(N), ')']);
legend('Calculated', 'Reference Values', 'Location', 'best');
saveas(gcf, filename11);


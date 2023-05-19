function [p] = pressure_plot(p, tx, th, N)

points = tx(1:end-1) + 0.5*th;
[X,Y] = meshgrid(points,points);

p = [0; p(1:N); 0; p(N+1:end-N); 0; p(end-N+1:end); 0];
p = reshape(p, N+2, N+2)';

p = p(2:end-1, 2:end-1);

% Y = flip(Y, 1);

a = 0.3;
b = 0.17;
c = 0.12;
d = 0.11;
e = 0.09;
f = 0.07;
g = 0.05;
h = 0.02;
i = 0;
j = -0.002;

contour(X, Y, p)%, 'LevelList',[a b c d e f g h i j])


end
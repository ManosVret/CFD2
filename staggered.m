function [X,Y] = staggered(tx)
x1 = tx(1:end-1)
x2 = tx(2:end)
avgX = (x1 + x2)/2
[X,Y] = meshgrid(avgX,avgX);
end
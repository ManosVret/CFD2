function [Ht02] = hodget02(h, n)
n_areas = (n+1)*(n+1);
areas = [];
for i=1:n+1
   for j = 1:n+1
        areas = [areas, h(i)*h(j)];
   end
end

Ht02 = sparse(diag(1./areas));

end
function [tE10] = maketE10(n)
edges = 2*(n*(1+n)); %Rows
points = (n+1)*(n+1); %Columns

%%Top half of tE10 matrix
diag1 = ones(1, edges/2 + n+1)*-1;
diag2 = ones(1, edges/2);

diag1 = diag(diag1);
diag2 = diag(diag2, n+1);

topblock = diag1+diag2;

topblock(end-n:end,:) = [];


%%Bottom half of tE10 matrix
diag1 = ones(1,n+1);
diag2 = -1*ones(1,n);
diag1 = diag(diag1);
diag2 = diag(diag2,1);

coreblock = diag1 + diag2;
coreblock(end,:) = [];

blklst = {};
for i = 1:n+1
    blklst{end+1} = coreblock;
end

bottomblock = blkdiag(blklst{:});


tE10 = sparse([topblock;bottomblock]);



end
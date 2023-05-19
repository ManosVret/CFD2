function [tE21, boundary, u_norm] = maketE21(n)
volumes = n*n+4*n;
edges = 2*n*(n+3);

%%Generate left side of the tE21 matrix (horizontal fluxes)
diag1 = ones(1,n+3)*-1;
diag2 = ones(1,n+2);
diag1 = diag(diag1);
diag2 = diag(diag2,1);
coreblock = diag1 + diag2;
coreblock(end,:) = [];

blocksize = size(coreblock);

zb_height = blocksize(1);
zb_width = blocksize(2);

coreblock = coreblock;

leftside = [coreblock, zeros(zb_height, zb_width*(n-1))];

for i = 2:n-1
    lzero = zeros(zb_height, zb_width*(i-1));
    rzero = zeros(zb_height, zb_width*(n-i));
    leftside = [leftside;lzero,coreblock,rzero];
end

leftside = [leftside;zeros(zb_height, zb_width*(n-1)),coreblock];


topblock = zeros(n,edges/2);

leftside = [topblock; leftside; topblock];

clear topblock
clear diag1
clear diag2

%%Generate right side of tE21 matrix (vertical fluxes)

diag1 = ones(1,n)*-1;
diag2 = ones(1,n);

diag1 = diag(diag1);
diag1 = [diag1, zeros(size(diag1))];
diag2 = diag(diag2, n);
diag2(end-n+1:end, :) = [];
coreblock = diag1 + diag2;
bufferrow = zeros(1,edges/2);
topzblock = zeros(n,edges/2-2*n);

rightside = [coreblock,topzblock;bufferrow;zeros(n,n),coreblock,zeros(n,edges/2 - 3*n);bufferrow];


for i = 3:n
    lzblock = zeros(n,n*(i-1));
    rzblock = zeros(n,edges/2 - n*(i+1));
    rightside = [rightside;bufferrow;lzblock,coreblock,rzblock;bufferrow];
end
rightside = [rightside;bufferrow;zeros(n,edges/2 - 3*n),coreblock,zeros(n,n);bufferrow;topzblock,coreblock];


btE21 = [leftside, rightside];  %Full matrix tE21 including boundaries

%%Removing boundary columns
tE21 = [];%Final tE21 matrix
bdy = [];%Matrix with boundary points


bdy = [bdy, btE21(:,1)];
idx = 2;

%Horizontal fluxes
for i=1:n-1
    for j=1:n+1
        tE21 = [tE21, btE21(:,idx)];
        idx = idx+1;
    end
    bdy = [bdy, btE21(:,idx)];
    idx = idx + 1;
    bdy = [bdy, btE21(:,idx)];
    idx = idx + 1;
end


for j=1:n+1
    tE21 = [tE21, btE21(:,idx)];
    idx = idx+1;
end

%vertical fluxes
for i = 1:n+1;
    bdy = [bdy, btE21(:,idx)];
    idx = idx + 1;
end

for i = 1:(n*(n+1))
    tE21 = [tE21, btE21(:,idx)];
    idx = idx+1;
end

for i = 1:n
    bdy = [bdy, btE21(:,idx)];
    idx = idx + 1;
end


tE21 = sparse(tE21);
boundary = sparse(bdy);
u_norm = boundary*zeros(4*n,1);
end
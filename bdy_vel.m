% bdy_vel(0, -1, 0, 0, [0.25,0.5,0.25, 0.35], 3)
function [u_pres] = bdy_vel(v_bot, v_top, v_left, v_right, h, n)

surf = (n+1)*(n+1);
bdy_edg = 4*(n+1);


% Boundary velocity vector, scaled with length of 
V_bdy = [repmat(v_bot, n+1, 1); repmat(v_top, n+1, 1); repmat(v_left, n+1, 1); repmat(v_right, n+1, 1)];

V_bdy = sparse(V_bdy.*(repmat(h, 1, 4)'));

%Calculation of contribution of circulation from these velocities on the
%dual grid surfaces

%left side of matrix
coreblock = diag(ones(1,n+1));
adj_zblock = zeros(n+1,n+1);

midblock = zeros(surf-2*(n+1),2*(n+1));

leftmat = [coreblock, adj_zblock; midblock; adj_zblock, -1.*coreblock];


%rights side of matrix
coreblock = zeros(n+1,2);
coreblock(1,1) =  1;
coreblock(2,2) = -1;

blklst = {};
for i = 1:n
    blklst{end+1} = coreblock;
end

diagblock = blkdiag(blklst{:});
zblock = zeros((n+1)*n,1);
midblock = [zblock, diagblock, zblock];

bottomblock = zeros(1, bdy_edg/2);
bottomblock(end) = 1;

topblock = zeros(n, bdy_edg/2);
topblock(1,1) = -1;

rightmat = [topblock; midblock; bottomblock];

%synthesise full matrix and find vorticity from prescribed tangential
% velocities
bE21 = sparse([leftmat, rightmat]);

u_pres = bE21*V_bdy;

end
function output = Skeleton_NS_solver(N, Re, tol, K, CFL)
    % The system that you need to solve will be singular. Matlab gives you a
    % warning at each time step. To switch off this warning, remove the comment
    % in the next line

    warning on

    % 00D#MMXXI#

    % This file contains the skeleton of the program which will solve the lid
    % driven cavity problem on a unit square. The parts that have to be
    % supplemented are described in the assignment.
    %
    % The pieces that need to be filled in are indicated
    %

    %
    % When running the code, determine a suitable time step. A too small time
    % step will make the calculation very long, while for a too large time step
    % the solution will blow up due to numerical instability.
    %

    Delta = 1/N;            % uniform spacing to be used in the mapping to compute tx

    % Determine a suitable time step and stopping criterion, tol

    % wall velocities
    U_wall_top = -1;
    U_wall_bot = 0;
    U_wall_left = 0;
    U_wall_right = 0;
    V_wall_top = 0;
    V_wall_bot = 0;
    V_wall_left = 0;
    V_wall_right = 0;

    %
    %   Generation of a non-uniform mesh
    %

    %
    %   tx are the coordinates of the nodal points on the outer-oriented
    %   primal mesh
    %

    tx = zeros(1,N+1);
    for i=1:N+1
        xi = (i-1)*Delta;
        tx(i) = 0.5*(1. - cos(pi*xi));
    end

    % Mesh width on the outer-oriented primal mesh
    th = zeros(N,1);
    th = tx(2:N+1) - tx(1:N);

    %
    %  x are the coordinates of the nodal points on the dual inner-orineted
    %  mesh (including endpoints 0 and 1)
    %  h contains the edge lengths on the inner-oriented dual mesh
    %
    x = 0.5*(tx(1:N) + tx(2:N+1));
    x = [0 x 1];

    h = zeros(N+1,1);
    h = x(2:N+2) - x(1:N+1);

    % Stable time step. Note that this is a conservative estimate, so it is
    % possible to increase this time step somewhat.

%     if nargin < 4
%         dt = min(min(h),0.5*Re*min(h)^2);
%     end

    
    dt = K * min(min(h),0.5*Re*min(h)^2);
    
    %
    %   Initial condition u=v=0
    %
    %   Both u and v will be stored in one big vector called 'u'
    %
    %   The vector u only contains the true unknowns, not the velocities that
    %   are prescribed by the boundary conditions
    %
    %   The vector u contains the *inner-oriented* circulations as unknowns

    u = zeros(2*N*(N+1),1);
    
    hmesh = zeros(size(u));
    for i=1:N
        hmesh([(i-1)*(N+1)+1:(i-1)*(N+1)+1+N]) = h;
    end
    for i=1:N+1
        hmesh([N*(N+1)+1+(i-1)*N:N*(N+1)+1+(i-1)*N+N-1]) = h(i);
    end
    % Set up the Incindence matrix 'tE21' which connects the fluxes to the
    % volumes. Use the orientation described in the assignment.

    %
    %  Inserting boundary conditions for normal velocity components and store
    %  this part in the vector u_norm, see assignemnt.
    %
    [tE21, boundary, u_norm] = maketE21(N);

    %  Set up the sparse, outer-oriented incidence matrix tE10.
    [tE10] = maketE10(N);

    %  Set up the sparse, inner-oriented  incidence matrix E10
    E10 = -tE21.';

    %  Set up the (extended) sparse, inner-oriented incidence matrix E21
    E21 = tE10.';

    %  Split off the prescribed tangential velocity and store this in
    %  the vector u_pres
    u_pres = bdy_vel(U_wall_bot, U_wall_top, V_wall_left, V_wall_right, h, N);

    %  Set up the Hodge matrices Ht11 and H1t1
    [Ht11, H1t1] = hodges11(h, th, N);

    %  Set up the Hodge matrix Ht02
    [Ht02] = hodget02(h, N);
    H2t0 = spdiags(1 ./ spdiags(Ht02, 0), 0, size(Ht02, 1), size(Ht02, 2)); % for int vort
    %
    % The prescribed velocities will play a role in the momentum equation
    %
    u_pres_tvort = Ht02 * u_pres; % U_pres to outer oriented 0 form representing contribution of boundary conditions to point-wise vorticity
    u_pres = H1t1 * E21' * Ht02 * u_pres; % U_pres to inner oriented 1 forms

    % Now all matrices are set up and the time stepping can start. 'iter' will
    % record the number of time steps. This allows you to give output after a
    % preselected number of time steps.
    %
    % 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
    % small, steady state has been reached. Determine a suitable value for
    % 'tol'
    %

    convective = zeros(2*N*(N+1),1);
    ux_xi = zeros((N+1)*(N+1),1);
    uy_xi = zeros((N+1)*(N+1),1);

    % Set up the matrix for the Poisson equation
    A = -tE21 * Ht11 * tE21';

    % Perform an LU-decomposition for the pressure matrix A
    [L, U] = lu(A);

    % Abbreviation for some matrix products which are constant in the time loop
    VLaplace = H1t1 * E21' * Ht02 * E21;
    DIV = tE21 * Ht11;
    
    diff = 1;
    iter = 1;
    diff_list = [];

    while diff > tol
    %while iter < 10
        
%         % Implement CFL timestep after first 10 iters
%         if iter > 10 && nargin > 3
%             dt = CFL * min(h) / max(abs(umesh));
%         end

        % Vector xi is obtained. It corresponds with the point-wise vorticity
        % at each cell

        % Vectors ux_xi and uy_xi correspond with the multiplications of
        % xi with the horizontal and vertical velocity components at each cell.
        % Only the cells required for vector convective are calculated. The
        % ordering of each vector with respect to the ordering of cells in the
        % grid is different (left to right for ux_xi and bottom to top for
        % uy_xi)

        xi = Ht02 * E21 * u + u_pres_tvort;
        xi_other = H2t0 * xi;
        
        for i = 1:N+1
            for j = 1:N+1
                k = j + (i-1)*(N+1);
                if j==1
                    ux_xi(k) = U_wall_bot * xi(i+(j-1)*(N+1));
                    uy_xi(k) = V_wall_left * xi(j+(i-1)*(N+1));
                elseif j==N+1
                    ux_xi(k) = U_wall_top * xi(i+(j-1)*(N+1));
                    uy_xi(k) = V_wall_right * xi(j+(i-1)*(N+1));
                else
                    ux_xi(k) = (u(i+(j-1)*(N+1)) + u(i+(j-2)*(N+1))) * xi(i+(j-1)*(N+1)) / (2. * h(i));
                    uy_xi(k) = (u(N*(N+1)+j+(i-1)*N) + u(N*(N+1)+j-1+(i-1)*N)) * xi(j+(i-1)*(N+1)) / (2. * h(i));
                end
            end
        end

        for i = 1:N
            for j = 1:N+1
                convective(j+(i-1)*(N+1)) = -(uy_xi(j+(i-1)*(N+1)) + uy_xi(j+i*(N+1))) * h(j) / 2.;
                convective(N*(N+1)+i+(j-1)*N) = (ux_xi(j+(i-1)*(N+1)) + ux_xi(j+i*(N+1))) * h(j) / 2.;
            end
        end

        % Set up the right hand side for the Poisson equation for the pressure
        rhs_Poisson = DIV * (u/dt - convective - VLaplace*u/Re - u_pres/Re) + u_norm/dt;

        % Solve for the new pressure
        temp = L \ rhs_Poisson;
        p = U \ temp;

        % Store the velocity from the previous time step in the vector u_old
        uold = u;

        % Update the velocity field
        u = u - dt * (convective - tE21'*p + VLaplace*u/Re + u_pres/Re);

        %
        % Every other 1000 iterations check whether you approach steady state
        % and check whether you satisfy conservation of mass. The largest
        % rate at which mass is destroyed or created is denoted by 'maxdiv'.
        % This number should be very small, in the order of machine precision.
        %
        umesh = u./hmesh;
        
        if mod(iter,1000) == 0
            maxdiv = max(DIV * u + u_norm);
            diff = max(abs(u - uold)) / dt
        end
        diff_list = [diff_list, diff];
        iter = iter + 1;
    end
    
    %%% Generate and format outputs %%%
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
    
    yref = flip([1.00000 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5000 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547]);
    uref = flip([-1.0000000 -0.6644227 -0.5808359 -0.5169277 -0.4723329 -0.3372212 -0.1886747 -0.0570178 0.0620561 0.1081999 0.2803696 0.3885691 0.3004561 0.2228955 0.2023300 0.1812881]);
    xref = [ 0.0312 0.0391 0.0469 0.0547 0.0937 0.1406 0.1953 0.5000 0.7656 0.7734 0.8437 0.9062 0.9219 0.9297 0.9375 1.0000];
    vref = [ -0.2279225 -0.2936869 -0.3553213 -0.4103754 -0.5264392 -0.4264545 -0.3202137 0.0257995 0.3253592 0.3339924 0.3769189 0.3330442 0.3099097 0.2962703 0.2807056 0.0000000];

    interp_uxstag = spline(cumsum(th), uxstag(:, floor(N/2))', yref);
    a = interp_uxstag - uref;
    error = rms(interp_uxstag - uref);

    % Specify the desired outputs in a structure
    output.u = u;
    output.p = p;
    output.xi = xi;
    output.diff_list = diff_list;
    output.iter = iter;
    output.h = h;
    output.th = th;
    output.error = error;
    output.interp_uxstag = interp_uxstag;
    output.ustag = ustag;
    output.A = A;
    output.h = h;
    output.th = th;
    output.x = x;
    output.xiother = xi_other;
end

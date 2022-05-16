clear
N = 3;

% Size init
Y = (N+2)^2 - 4;
X = 2*N^2 + 6*N;
A = zeros( Y , X );

for i = 1:N
    % Cells with vertical fluxes
    A(i,i) = -1;
    A(i,i+N) = 1;
    
    A(Y+1-i,X+1-i) = 1;
    A(Y+1-i,X+1-i-N) = -1;
    
    % Cells with horizontal fluxes
    A((N+1) + (N+2)*(i-1) , (2*N+1) + (2*N+3)*(i-1)) = -1;
    A((N+1) + (N+2)*(i-1) , (2*N+1) + (2*N+3)*(i-1) + 1) = 1;
    
    A((N+1) + (N+2)*(i-1) +N+1, (2*N+1) + (2*N+3)*(i-1) +N+1) = -1;
    A((N+1) + (N+2)*(i-1) +N+1, (2*N+1) + (2*N+3)*(i-1) + 1 +N+1) = 1; 
    
    % Cells with all 4 fluxes
    posleft = (N+2)*i;
    for j = 1:N
        pos = posleft + j-1;
        bottom = (N+1) + (2*N+3)*(i-1) +j-1;
        A(pos, bottom) = -1;
        A(pos, bottom + N+1) = -1;
        A(pos, bottom + N+2) = 1;
        A(pos, bottom + 2*N + 3) = 1;
    end
end




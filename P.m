function P = P(p,N)
I = 1;
J = 1;
P=[];
count = 1;
%for loopt over the columns
for i = 1:(N+2)
    %first columns is only boundarypoints
    if i == 1
        for j = 1:N
        count = count +1;
        end

    %middle columns consist of inner points and two boundary at each end
    elseif i == N+2
        for j = 1:N
            count = count +1;
        end

    else
        for j = 1:(N+2)
            if j == 1
                count = count +1 ;
            elseif j== N+2
                count = count +1 ;
            else 
                P(j-1,i-1) = p(count);
                count = count +1;
            end
        end
    end
    

end
P = P - P(round((N+1)/2),round((N+1)/2));
end


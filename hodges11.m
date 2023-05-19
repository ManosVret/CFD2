function [Ht11, H1t1] = hodges11(h,th, n)

    te = repmat(th, 1, n+1);
    e= repmat(h,1,n);
    tend = [];
    
    for i=1:n+1
        e = [e, repmat(h(i), 1, n)];
    end
    
    for i=1:n
        tend = [tend, repmat(th(i), 1, n+1)];
    end
    te = [tend, te];

    H1t1 = sparse(diag(e./te));
    Ht11 = sparse(diag(te./e));
end
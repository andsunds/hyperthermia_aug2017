function [ A ] = F_coef_mtrx(R, X, Y )

%tri = @(m) m*(m+1)/2;
L=size(X,1);
A=zeros(L);

for m = 1:L
    for n = 1:m
        A(m,n) = F_coef(R, X(m,n), Y(m,n), 0 );
    end
end

A=A+A.'-diag(diag(A));



end


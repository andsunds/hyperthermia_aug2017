function [ D, V ] = ring2D( L, c, q, rot_order, n_eigs )
% Calculates the largest eigenvalues for circle rings, using the Bessel method, in 2 dimensions.
%
% L:         size of the discretization (LxL matrix)
% c:         scaling factor between the frequency and space disc
% rot_order: the rotational order of the eigenfunction (0 for largest eigenvalue)
% n_eigs:    number of eigenvalues tobe calculated
%
% D: eigenvalue(s), normalized to be the eigevalue of the original problem
% V: corresponding eigenfunction(s)



%           Rs=1
v=c*linspace(q,1,L)'; %=k*Rs
V2=repmat(v',L,1);

Ks=coef_mtrx( v, rot_order ).*V2;




if nargout >1
    [V,D]=eigs(Ks,n_eigs,'LM');  %eigenfunction(s) and eigenvalue(s)
    V=V./repmat(sqrt(u),1,n_eigs); %Normalizing the V to the eigenfunctions of the original problem
else
    D    =eigs(Ks,n_eigs,'LM');  %only eigenvalue(s)
end
D=D*c*(1-q)/L; %/Rs^2, but we have difined Rs:=1.


end


function [ A ] = coef_mtrx( v, rot_order )

L=length(v); %We get the discretization size from x.

tri=@(i) i.*(i+1)/2; %the i'th triangle number
%An indexed lower triangle of a matrix together with row and col numbers.
%   m:  _  
%   1  |1|_  :j
%   2  |2 3|_ 
%   3  |4 5 6|__
%   4  |7_8_9_10|
%    n: 1 2 3 4
%Given an index j, the row (m) and col (n) can be calculated as below. To
%calculate m, all we have to do is invert the triangle number formula
J = 1:tri(L); %all indeces
m = ceil(-.5 + sqrt(.25 +2*J));
n = J -tri(m-1);

%For each j, we need the values x(m(j)) and y(n(j)). To optimize the
%parallellization, new (and longer) vectors X and Y are created precisely
%so that X(j) = x(m(j)) and Y(j) = y(n(j)).
v1=v(m); 
v2=v(n); %x(n) is sufficient since x and y are discretized identically

a=zeros(1,tri(L)); %init
%Parallellized for loop over all combinations of x and y.
for j=J
    if m(j)==n(j)
        a(j)=.5*((besselj(rot_order,v1(j))).^2-...
                 besselj(rot_order-1,v1(j)).*besselj(rot_order+1,v1(j)));
    else
        a(j)=(v2(j)*besselj(rot_order,v1(j))*besselj(rot_order-1,v2(j))-...
              v1(j)*besselj(rot_order-1,v1(j))*besselj(rot_order,v2(j)))/...
              (v1(j)^2-v2(j)^2);
    end
end

%Transfering the parallellized (1D) results to the matrix A
A=zeros(L);
for j=J
    A(m(j),n(j)) = a(j);
end
%A is symmetric, which is why we only calculate half of A and then flip it
%over on itself; if just A+A', we get double the value on the diagonal
%which is why -diag(diag(A)) is needed. 
A=A+A.'-diag(diag(A));
% note: important to use .'. If rot_order =/= 0, A is not hermetian, only
% symmetric.

end
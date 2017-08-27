function [ A ] = PAR_F_coef_mtrx( R, k, rot_order )
%   PARALLEL version of F_coef_mtrx(R, x,y )
%
%   Returns all the Fourier coefficients for each pair of k and kappa.
%
%In the problem of finding the optimal frequency distribution in a circle
%ring of finite thickness for maximum energy in a circle in real space, we
%need to find the eigenvalue of
%   \lambda b_m(k) 
%   = R^2\int_{q\Gamma}^{\Gamma} \rd\kappa a_m(k,\kappa)b_m(\kappa).
%This function calculates a_0(k,\kappa). We only need to calculate m=0
%since that's the optimal case. 
%
% R:         size of the disc in real space (corresponds to c)
% x:         vector of discretized k values
% rot_order: rotational order (0 for largest eigenvalue)

L=length(k); %We get the discretization size from x.

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
k1=k(m); 
k2=k(n); %x(n) is sufficient since x and y are discretized identically

a=zeros(1,tri(L)); %init
%Parallellized for loop over all combinations of x and y.
parfor j=J
    %Calculating the Fourier coefficients.
    a(j) = F_coef(R, k1(j), k2(j), rot_order );
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


function [ a ] = F_coef( R, k1, k2, rot_order )
%   Numerically calculates the Fourier coefficients a_m(k1, k2).
%This is used in each loop iteration above.

% We need to find the coefficients
%   a_m(k1,k2) = 1/(4\pi) 
%                \int_0^{2\pi} [J_0(z)+J_2(z)]\ee^{-\ii m\theta} \id\theta,
% where z = R \sqrt{k1^2+k2^2 - 2k1k2 \cos(t)}.

z = @(t) R*sqrt(k1^2 + k2^2 - 2*k1*k2*cos(t));

a = 1/(4*pi)*integral(@(t) cos(rot_order*t).*...
                           (besselj(0,z(t)) + besselj(2,z(t))), 0, 2*pi);
%The only rotational dependence must be of cos form, since the kernel is
%even in t.
                       
end
function [ A ] = PAR_F_coef_mtrx( R, x, y )
%   PARALLEL version of F_coef_mtrx(R, x,y )
%   Returns all the Fourier coefficients for each pair of k and kappa.
%In the problem of finding the optimal frequency distribution in a circle
%ring of finite thickness for maximum energy in a circle in real space, we
%need to find the eigenvalue of
%   \lambda b_m(k) 
%   = R^2\int_{q\Gamma}^{\Gamma} \rd\kappa a_m(k,\kappa)b_m(\kappa).
%This function calculates a_0(k,\kappa). We only need to calculate m=0
%since that's the optimal case. 
%


L=size(x,1);
A=zeros(L);

%{
for m = 1:L
    for n = 1:m
        A(m,n) = F_coef(R, x(m,n), y(m,n), 0 );
    end
end
%}


tri=@(i) i.*(i+1)/2;
a=zeros(1,tri(L));
J = 1:tri(L);
m = ceil(-.5 + sqrt(.25 +2*J));
n = J -tri(m-1);

X=x(m);
Y=y(n);
parfor j=J
    a(j) = F_coef(R, X(j), Y(j), 0 );
end
for j=J
    A(m(j),n(j)) = a(j);
end


%A is symmetric, which is why we only calculate half of A and then flips it
%over on itself.
A=A+A.'-diag(diag(A));

end

function [ a ] = F_coef( R, k1, k2, m )
%   Numerically calculates the Fourier coefficients a_m(k1, k2).

% We need to find the coefficients
%   a_m(k1,k2) = 1/(4\pi) 
%                \int_0^{2\pi} [J_0(z)+J_2(z)]\ee^{-\ii m\theta} \id\theta,
% where z = R \sqrt{k1^2+k2^2 - 2k1k2 \cos(t)}.

z = @(t) R*sqrt(k1^2 + k2^2 - 2*k1*k2*cos(t));

a = 1/(4*pi)*integral(@(t) exp(-1i*m*t).*...
                           (besselj(0,z(t)) + besselj(2,z(t))), 0, 2*pi);

end
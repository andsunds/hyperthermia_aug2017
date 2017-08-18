function [ a ] = F_coef( R, k1, k2, m )
%   Numerically calculates the Fourier coefficients a_m(k1, k2).
%In the problem of finding the optimal frequency distribution in a circle
%ring of finite thickness for maximum energy in a circle in real space, we
%need to find the eigenvalue of
%   \lambda b_m(k) 
%   = R^2\int_{q\Gamma}^{\Gamma} \rd\kappa a_m(k,\kappa)b_m(\kappa).
%This function calculates a_m(k,\kappa).
%
% We need to find the coefficients
%   a_m(k1,k2) = 1/(4\pi) 
%                \int_0^{2\pi} [J_0(z)+J_2(z)]\ee^{-\ii m\theta} \id\theta,
% where z = R \sqrt{k1^2+k2^2 - 2k1k2 \cos(t)}.

z = @(t) R*sqrt(k1^2 + k2^2 - 2*k1*k2*cos(t));

a = 1/(4*pi)*integral(@(t) exp(-1i*m*t).*...
                           (besselj(0,z(t)) + besselj(2,z(t))), 0, 2*pi);

end


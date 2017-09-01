function [ D, V ] = ring2D_sym( L, c, q, rot_order, n_eigs )
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
Ks=calc_sym_mtrx( @(a,b) coef(a,b,rot_order), v);

if nargout >1
    [V,D]=eigs(Ks,n_eigs,'LM');  %eigenfunction(s) and eigenvalue(s)
else
    D    =eigs(Ks,n_eigs,'LM');  %only eigenvalue(s)
end
D=D*c*(1-q)/L; %/Rs^2, but we have difined Rs:=1.


end

function [a] = coef(v1,v2,rot_order)
%The coefficient function that calculates the coefficients of the symmetric
%matrix
if v1==v2
    a=.5*((besselj(rot_order,v1)).^2-...
          besselj(rot_order-1,v1).*besselj(rot_order+1,v1));
else
    a=(v2*besselj(rot_order,v1)*besselj(rot_order-1,v2)-...
       v1*besselj(rot_order-1,v1)*besselj(rot_order,v2))/...
       (v1^2-v2^2);
end
a=a*sqrt(v1*v2);
end

function [ D, V ] = ring2D_AB( L, c, q, headsize, rot_order, n_eigs )
% Calculates the largest eigenvalues for discs in 2 dimensions.
%
%Modified to get a SYMMETRIC kernel according to Slepian IV (1964) eqn. 
%(19) and (20).
%Saves around a factor 10 (!) in runtime compared to non-symmetrized.
% L:         size of the discretization (LxL matrix)
% c:         scaling factor between the frequency and space disc
% rot_order: the rotational order of the eigenfunction (0 for largest eigenvalue)
% n_eigs:    number of eigenvalues tobe calculated
%
% D: eigenvalue(s), normalized to be the eigevalue of the original problem
% V: corresponding eigenfunction(s)



%           Rs=1
vt=c*linspace(q,1,L)'; %v=kRs
vh=c*headsize*linspace(q,1,L)'; %v=kRs

%Ks=besselj(rot_order,c* XY ).*sqrt(c*XY); %sym. kernel according to (20)
H=calc_sym_mtrx( @(a,b) coef(a,b,rot_order), vh);
T=calc_sym_mtrx( @(a,b) coef(a,b,rot_order), vt);






D=eigs(T,H,n_eigs);  %only eigenvalue(s)

[~,I]=max(abs(D));
%Think about the normalization!
D=D(I)*c*(1-q)/L*headsize^2; %/Rs^2, but we have difined Rs:=1.


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


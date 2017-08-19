function [ DD, V ] = disc2D_sym( L, c, rot_order, n_eigs  )
% Calculates the largest eigenvalue for discs.
%Modified to get a SYMMETRIC kernel according to Slepian IV (1964) eqn. 
%(19) and (20).
%Saves around a factor 10 (!) in runtime.
% L:         size of the discretization (LxL matrix)
% c:         scaling factor between the frequency and space disc
% rot_order: the rotational order of the eigenfunction (0 for largest eigenvalue)
% n_eigs:    number of eigenvalues tobe calculated
%
% DD: eigenvalue, normalized to be the eigevalue of the original problem
% V:  eigenfunction


eps=1/L; %can't have x or y being 0
% For some reason 1/L seems to be the sweet spot. If eps < 1/L the
% eigenvalues get too small (comparison with Sleipan), and if eps>1/L the 
% eigenvalues become too large (greater than 1, which is theoretically
% impossible).

x=linspace(eps,1,L)';   %the free time variable of kernel
%y=linspace(eps,1,L)';  %integrated freq variable in the kernel
% We don't need a separet vecor y since both x and y live on the interval 
% 0 to 1, and since the kernel is symmetric.
XY=x*x'; %This creates a matrix whose element (i,j) is x(i)*y(j)

Ks=besselj(rot_order,c* XY ).*sqrt(c*YX); %sym. kernel according to (20)

[V,D]=eigs(Ks,n_eigs,'LM');  %eigenfunctions and eigenvalues
D=diag(D)/L;                 %vector form, and normalizing

d=2; %2 dimensional disc
DD=(c/2/pi)^d.*abs(2*pi*D).^2/c; %have to change normalization to get the wanted eigenvalues

end


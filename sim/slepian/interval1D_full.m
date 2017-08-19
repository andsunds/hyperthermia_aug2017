function [ D, V, x ] = interval1D_full( L, c, n_eigs )
% Calculates the largest eigenvalues for full intervals in 1 dimension.
%
%The integral eigenvalue equation
%   \lambda \phi(x) = \int K(x,y) \phi(y) dy
%is discretized via
%   \int dy  \to  \sum_{n} \frac{1}{N},
%which means that that eigenvalue equation becomes
%   \lambda \phi_m = \frac{1}{N} \sum_n K_{m,n} \phi_n.
%This is just a regular matrix eigenvalue equation. This is what we're
%implementing here. 
%Note that the $\lambda$ we're seeking is a factor % $T/N$ times the 
%eigenvalue of the matrix $K_{m,n}$.
%
% L:         size of the discretization (LxL matrix)
% c:         scaling factor between the frequency and space interval
%            c=2*W*T to agree with the notation of Slepian I (1961)
% n_eigs:    number of eigenvalues tobe calculated
%
% D: normalized eigenvalues
% V: corresponding eigenfunctions




T=1;                                %timespan
x=linspace(-T/2,T/2,L)';            %the free time variable of kernel
%y=linspace(-T/2,T/2,L)';           %integrated variable in the kernel
% no separate y, since x and y are identical in their discretizations.
xy=(repmat(x,1,L)-repmat(x',L,1));  %combine x and y into a matrix of x-y

W=2*c/T;                %frequency
Ks=sin(W*xy)./xy/pi;    %kernel: \frac{\sin(\Omega(s-t))}{\pi(s-t)}
%In the expression for the kernel above, we have some elements whic are NaN
%due to "0/0". If we look at the expression for Ks, we see that the NaN
%elements should be:
%   sin(W*(s-t))/(pi*(s-t)) --> W/pi,
%when (s-t)-->0.
Ks(isnan(Ks))=W/pi;

if nargout >1
    [V,D]=eigs(Ks,n_eigs,'LM');  %eigenfunction(s) and eigenvalue(s)
else
    D    =eigs(Ks,n_eigs,'LM');  %only eigenvalue(s)
end

D=diag(D)*T/L;    %vector form, and normalizing to the wanted eigenvalues


end


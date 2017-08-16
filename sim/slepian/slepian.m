%% D=1
% The integral eigenvalue equation
%   \lambda \phi(x) = \int K(x,y) \phi(y) dy
% is discretized via
%   \int dy  \to  \sum_{n} \frac{1}{N},
% which means that that eigenvalue equation becomes
%   \lambda \phi_m = \frac{1}{N} \sum_n K_{m,n} \phi_n.
% This is just a regular matrix eigenvalue equation. This is what we're
% implementing here. 
% Note that the $\lambda$ we're seeking is a factor % $1/N$ times the 
% eigenvalue of $K_{m,n}$.
clc;clear all; clf


L=5;                              %prec of int approx
T=1;                                %timespan
x=linspace(-T/2,T/2,L)';            %the free time variable of kernel
y=linspace(-T/2,T/2,L)';            %integrated variable in the kernel
xy=(repmat(x,1,L)-repmat(y',L,1));  %combine x and y into a matrix of x.-y

W=8;                    %freq
Ks=sin(W*xy)./xy/pi;    %kernel
Ks(isnan(Ks))=W/pi;     %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0


[V,D]=eig(Ks);  %eigenfunctions and eigenvalues
D=diag(D)/L;    %vector form, and normalizing to the wanted eigenvalues

i=find(D>1e-4);                 %finds the largest eigenvlaues
fprintf('D = %1.4d,\t\n',D(i)); %displays those eigenvalues
plot(V(:,i));                   %plots correspoding eigenfunctions





%% D=2, p=0
clc;clear all; clf


L=500;                  %prec of int approx
eps=.1/L;               %can't have y=0
x=linspace(eps,1,L)';   %the free time variable of kernel
y=linspace(eps,1,L)';   %integrated freq variable in the kernel
[Y,X]=meshgrid(x,y);    %combine x and y into mesh (note order: [Y,X])

c=1;                        % c = \Omega T
N=0;
Ks=besselj(N,c*X.*Y).*Y;    %kernel
%Ks(isnan(Ks))=W/pi;        %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0
%Ks

[V,D]=eig(Ks);  %eigenfunctions and eigenvalues
D=diag(D)/L;    %vector form, and normalizing


d=2;                            %dimensions
DD=(c/2/pi)^d.*abs(2*pi*D).^2;  %have to change normalization to get the wanted eigenvalues
i=find(DD>1e-4);                %finds the largest eigenvlaues

fprintf('DD = %1.4d,\t\n',sort(DD(i),1,'descend')); %displays those eigenvalues
plot(V(:,i));                                       %plots correspoding eigenfunctions














































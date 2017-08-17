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


L=200;                               %prec of int approx
T=1;                                %timespan
x=linspace(-T/2,T/2,L)';            %the free time variable of kernel
y=linspace(-T/2,T/2,L)';            %integrated variable in the kernel
xy=(repmat(x,1,L)-repmat(y',L,1));  %combine x and y into a matrix of x.-y

W=8;                    %freq
Ks=sin(W*xy)./xy/pi;    %kernel
Ks(isnan(Ks))=W/pi;     %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0

tic
[V,D]=eig(Ks);  %eigenfunctions and eigenvalues
D=diag(D)*T/L;    %vector form, and normalizing to the wanted eigenvalues
toc

i=find(D>1e-4);                 %finds the largest eigenvlaues
fprintf('D = %1.4d,\t\n',sort(D(i), 'descend')); %displays those eigenvalues
plot(x,V(:,i));                   %plots correspoding eigenfunctions


%% D=1, Simpson's method
clc;clear all; clf
% ¡Significantly slower (factor 10), due to Kmat no longer being symmetric!

L=2*4;                              %prec of int approx
T=1;                                %timespan
x=linspace(-T/2,T/2,L)';            %the free time variable of kernel
y=linspace(-T/2,T/2,L)';            %integrated variable in the kernel
xy=(repmat(x,1,L)-repmat(y',L,1));  %combine x and y into a matrix of x.-y

W=8;                    %freq
Ks=sin(W*xy)./xy/pi;    %kernel
Ks(isnan(Ks))=W/pi;     %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0

tic
U=2*(1+(-1).^(1:L)) + (1-(-1).^(1:L));
U(1)=1; U(L)=1;

Kmat= T/(3*L) * repmat(U,L,1).*Ks;

[V,D]=eig(Kmat);  %eigenfunctions and eigenvalues
D=diag(D);    %vector form, and normalizing to the wanted eigenvalues
toc

i=find(D>1e-4);                 %finds the largest eigenvlaues
fprintf('D = %1.4d,\t\n',sort(D(i), 'descend')); %displays those eigenvalues
plot(x,V(:,i));                   %plots correspoding eigenfunctions



%% D=2, p=0
clc;clear all; clf


L=1000;                  %prec of int approx
eps=.1/L;               %can't have y=0
x=linspace(eps,1,L)';   %the free time variable of kernel
y=linspace(eps,1,L)';   %integrated freq variable in the kernel
[Y,X]=meshgrid(x,y);    %combine x and y into mesh (note order: [Y,X])

c=8;                        % c = \Omega T
N=0;
Ks=besselj(N,c*X.*Y).*Y;    %kernel
%Ks(isnan(Ks))=W/pi;        %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0
%Ks
tic
[V,D]=eig(Ks);  %eigenfunctions and eigenvalues
D=diag(D)/L;    %vector form, and normalizing
toc

d=2;                            %dimensions
DD=(c/2/pi)^d.*abs(2*pi*D).^2;  %have to change normalization to get the wanted eigenvalues
i=find(DD>1e-4);                %finds the largest eigenvlaues

fprintf('DD = %1.4d,\t\n',sort(DD(i),1,'descend')); %displays those eigenvalues
plot(x,V(:,i));                                       %plots correspoding eigenfunctions



%% D=2, p=0, even kernel
clc;clear all; clf
%Once again, modifying to get an EVEN kernel saves around a factor 10 (!)

L=1000;                  %prec of int approx
eps=.1/L;               %can't have y=0
x=linspace(eps,1,L)';   %the free time variable of kernel
y=linspace(eps,1,L)';   %integrated freq variable in the kernel
[Y,X]=meshgrid(x,y);    %combine x and y into mesh (note order: [Y,X])

c=8;                        % c = \Omega T
N=0;
Ks=besselj(N,c*X.*Y).*sqrt(c*Y.*X);    %kernel
%Ks(isnan(Ks))=W/pi;        %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0
%Ks
tic
[V,D]=eig(Ks);  %eigenfunctions and eigenvalues
D=diag(D)/L;    %vector form, and normalizing
toc

d=2;                            %dimensions
DD=(c/2/pi)^d.*abs(2*pi*D).^2/c;  %have to change normalization to get the wanted eigenvalues
i=find(DD>1e-4);                %finds the largest eigenvlaues

fprintf('DD = %1.4d,\t\n',sort(DD(i),1,'descend')); %displays those eigenvalues
plot(x, V(:,i)./repmat(sqrt(x),1,length(i)));                                       %plots correspoding eigenfunctions












































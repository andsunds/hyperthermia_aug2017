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
[V,D]=eigs(Ks,6,'LM');  %eigenfunctions and eigenvalues
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

L=10;                  %prec of int approx
eps=1/L;               %can't have y=0
% For some reason 1/L seems to be the sweet spot. If eps < 1/L the
% eigenvalues get too small (comparison with Sleipan), and if eps>1/L the 
% eigenvalues they become too large (greater than 1, which is theoretically
% impossible).
x=linspace(eps,1,L)';   %the free time variable of kernel
y=linspace(eps,1,L)';   %integrated freq variable in the kernel
[Y,X]=meshgrid(x,y);    %combine x and y into mesh (note order: [Y,X])

c=1;                        % c = \Omega T
N=0;
Ks=besselj(N,c*X.*Y).*sqrt(c*Y.*X);    %kernel
%Ks(isnan(Ks))=W/pi;        %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0
%Ks
tic
[V,D]=eigs(Ks,6,'LM',optimset('maxit',300));  %eigenfunctions and eigenvalues
D=diag(D)/L;    %vector form, and normalizing
toc

d=2;                            %dimensions
DD=(c/2/pi)^d.*abs(2*pi*D).^2/c;  %have to change normalization to get the wanted eigenvalues
i=find(DD>1e-4);                %finds the largest eigenvlaues

fprintf('DD = %1.7d,\t\n',sort(DD(i),1,'descend')); %displays those eigenvalues
plot(x, V(:,i)./repmat(sqrt(x),1,length(i)));                                       %plots correspoding eigenfunctions



%% D=2, p=0, even kernel, c loop
clc;clear all; clf
%Once again, modifying to get an EVEN kernel saves around a factor 10 (!)

L=1000;                  %prec of int approx
eps=1/L;               %can't have y=0
% For some reason 1/L seems to be the sweet spot. If eps < 1/L the
% eigenvalues get too small (comparison with Sleipan), and if eps>1/L the 
% eigenvalues they become too large (greater than 1, which is theoretically
% impossible).
x=linspace(eps,1,L)';   %the free time variable of kernel
y=linspace(eps,1,L)';   %integrated freq variable in the kernel
[Y,X]=meshgrid(x,y);    %combine x and y into mesh (note order: [Y,X])

n=100;
data=zeros(n,1);
C=linspace(0.1,10,n);                        % c = \Omega T
for j=1:n
    c=C(j);
    N=0;
    Ks=besselj(N,c*X.*Y).*sqrt(c*Y.*X);    %kernel

    [V,D]=eigs(Ks,6,'LM',optimset('maxit',1000));  %eigenfunctions and eigenvalues
    D=diag(D)/L;    %vector form, and normalizing

    d=2;                            %dimensions
    DD=(c/2/pi)^d.*abs(2*pi*D).^2/c;  %have to change normalization to get the wanted eigenvalues
    i=find(DD>1e-4);                %finds the largest eigenvlaues

    data(j)=max(DD(i));
end
plot(C,1-data)

%% plots and symptotics
% Here we do see some numerical error in the eigenvalue calculation since
% they are not following the asymptotic formula given by Slepian (IV).

clf;clc;clear all

C=linspace(0.1,16,100); 
data=load('lambda_c0.1-16.tsv','-ascii');
plot(C,1-data)
hold on

%I=(20:50).';
%A=[C(I)', ones(size(I))]\log(1-data(I));
%plot(C, exp(A(2)+A(1)*C))

plot(C,pi*8*C.*exp(-2*C))% asymptotic formula

set(gca,'xlim',[0,10],'yscale','log')




%% D=2, not the full disc
%In the problem of finding the optimal frequency distribution in a circle
%ring of finite thickness for maximum energy in a circle in real space, we
%need to find the largest eigenvalue of
%   \lambda b_m(k) 
%   = R^2\int_{q\Gamma}^{\Gamma} \rd\kappa a_m(k,\kappa)b_m(\kappa).

clc;clear

R=2;
Gamma=1;
q=.1;

L=30;

x=Gamma*linspace(q,1,L)';   
y=Gamma*linspace(q,1,L)';   
[Y,X]=meshgrid(x,y);    


A = F_coef_mtrx(R, X, Y );

D=eigs(A.*Y,1)*R^2*Gamma*(1-q)/L



%TODO:
%See if Simpson's method helps here
%then do a c loop (maybe also a q loop).























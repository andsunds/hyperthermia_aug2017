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


L=5000; %prec of int approx
c=2;
n_eigs = 6;   % # of eigenvalues to be calc'ed for this rot order
tic
D=interval1D_full( L, c, n_eigs );
toc
fprintf('D = %1.4d,\t\n',sort(D, 'descend')); %displays those eigenvalues
%plot(x,V(:,i));                   %plots correspoding eigenfunctions


%% D=1, Simpson's method   >>> !NO GOOD! <<<
%!!!Significantly slower (factor 10), due to Kmat no longer being symmetric!!!
% Also not as good results (compared to Slepian I (1961) Table I).
%{
clc;clear all; clf

L=2*400;                              %prec of int approx
T=1;                                %timespan
x=linspace(-T/2,T/2,L)';            %the free time variable of kernel
y=linspace(-T/2,T/2,L)';            %integrated variable in the kernel
xy=(repmat(x,1,L)-repmat(y',L,1));  %combine x and y into a matrix of x.-y

W=2;                    %freq
Ks=sin(W*xy)./xy/pi;    %kernel
Ks(isnan(Ks))=W/pi;     %sin(W*(x-y))/(pi*(x-y)) -> W/pi when (x-y)->0

tic
% Simpson's method implemeted here becomes, with N = 2M,
%   \int_{k_0}^{k_N} \rd{k} f(k) 
%   \approx \frac{\Delta{k}}{3} \sum_{j=0}^{M-1}[f_{2j}+4f_{2j+1}+f_{2j+2}]
%   = \Delta{k}/3 [f_0 + 4f_1 + 2f_2 + 4f_3 + 2f_4 + ... + 4f_{N-1} + f_N]
% The vector U represents this sum

U=2*(1+(-1).^(1:L)) + (1-(-1).^(1:L));
U(1)=1; U(L)=1;

Kmat= T/(3*L) * Ks*diag(U);

[V,D]=eig(Kmat);  %eigenfunctions and eigenvalues
D=diag(D);    %vector form, and normalizing to the wanted eigenvalues
toc

i=find(D>1e-4);                 %finds the largest eigenvlaues
fprintf('D = %1.4d,\t\n',sort(D(i), 'descend')); %displays those eigenvalues
plot(x,V(:,i));                   %plots correspoding eigenfunctions
%}


%% D=2, asymmetric kernel   >>> !NO GOOD! <<<
% The most straight forward way to find the eigenvalues for 2D discs, from
% Slepian IV (1964) eqn. (18).
%{
clc;clear all; clf

L=1000;                 %prec of int approx
eps=1/L;                %can't have y=0
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
plot(x,V(:,i));                                     %plots correspoding eigenfunctions

%}


%% D=2, symmetric kernel
clc;clear all; clf
% Modified to get a SYMMETRIC kernel according to Slepian IV (1964) eqn.
% (19) and (20).
% Saves around a factor 10 (!) in runtime.

L=100;                  %prec of int approx

c         = 10;  % c = \Omega T
rot_order = 0;   % rotational order of full eigenfunction (just use 0)
n_eigs    = 2;   % # of eigenvalues to be calc'ed for this rot order
tic
DD = disc2D_sym( L, c, rot_order, n_eigs );
toc


fprintf('DD = %1.7d,\t\n',sort(DD,1,'descend')); %displays those eigenvalues
%plot(x, V(:,i)./repmat(sqrt(x),1,length(i)));      %plots correspoding eigenfunctions



%% D=2, even kernel, c loop
clc;clear all; clf
% Same code as above but now running i a loop over different c values.

L=500; %prec of int approx

n=8;
data=zeros(n,1);
C=[.1 .5 1 1.5 2 3 5 10]; % c = \Omega T

rot_order = 0;   % rotational order of full eigenfunction (just use 0)
n_eigs    = 2;   % # of eigenvalues to be calc'ed for this rot order

tic
parfor j=1:n
    c=C(j);
    DD = disc2D_sym( L, c, rot_order, n_eigs );
    data(j)=max(DD);
end
toc
for j=1:n
    fprintf('c = %3.1f \t\tDD = %1.8e\n',C(j),data(j))
end


%% plots and asymptotics   (needs more love)
% Here we do see some numerical error in the eigenvalue calculation since
% they are not following the asymptotic formula given by Slepian IV (1964).


clf;clc;clear all

C=linspace(0.1,16,100); 
data=load('lambda_c0.1-16.tsv','-ascii');
plot(C,1-data)
hold on

%I=(20:50).';
%A=[C(I)', ones(size(I))]\log(1-data(I));
%plot(C, exp(A(2)+A(1)*C))

plot(C,pi*8*C.*exp(-2*C))%asymptotic formula eqn. (93) of Slepian IV (1964) 

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

L=3;

x=Gamma*linspace(q,1,L)';   
%y=Gamma*linspace(q,1,L)'; %non-parallellized
%[Y,X]=meshgrid(x,y);      %non-parallellized
Y=repmat(x',L,1);

tic
A = PAR_F_coef_mtrx(R, x);
%A = F_coef_mtrx(R, X, Y); %non-parallellized
toc

tic
D=eigs(A.*Y,1)*R^2*Gamma*(1-q)/L;
toc



%% D=2, not the full disc, c-loop 
%In the problem of finding the optimal frequency distribution in a circle
%ring of finite thickness for maximum energy in a circle in real space, we
%need to find the largest eigenvalue of
%   \lambda b_m(k) 
%   = R^2\int_{q\Gamma}^{\Gamma} \rd\kappa a_m(k,\kappa)b_m(\kappa).

clc;clear
L=1000;
q=1/L; %q determins the inner radius of the ring (1/L means the full disc)

n=8;
C=[.1 .5 1 1.5 2 3 5 10];

R=1;
A=cell(n,1);     %init
Gamma=cell(n,1); %init
x=cell(n,1);     %init


tic
for i=1:n
    c=C(i);
    Gamma{i}=c/R;
    x{i}=Gamma{i}*linspace(q,1,L)';   

    A{i} = PAR_F_coef_mtrx(R, x{i}, 0);
end
toc

tic
data=zeros(n,1); %init
parfor j=1:n
    Y=repmat(x{j}',L,1);
    data(j)=eigs(A{j}.*Y,1)*R^2*Gamma{j}*(1-q)/L;
end
toc

for j=1:n
    fprintf('c = %3.1f \t\tDD = %1.8e\n',C(j),data(j))
end




















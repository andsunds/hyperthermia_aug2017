%% D=2, even kernel, c loop
clc;clear all;
% Same code as above but now running i a loop over different c values.

L=10000; %prec of int approx

n=101;
data=zeros(n,2);
%C=[.1 .5 1 1.5 2 3 5 10]; % c = \Omega T
C=linspace(0,10,n);

rot_order = 0;   % rotational order of full eigenfunction (just use 0)
n_eigs    = 1;   % # of eigenvalues to be calc'ed for this rot order

tic
parfor j=1:n
    %fprintf('n = %3d\n',j)
    c=C(j);
    D1 = interval1D_full( L, c, n_eigs );
    D2 = disc2D_sym( L, c, rot_order, n_eigs );
    data(j,:)=[D1,D2];
end
time=toc;
hr  = floor(time/3600); time=time-hr*3600;
min = floor(time/60); time=time-min*60;
sec=floor(time); 
secdec=floor(1000*(time-sec));

file=fopen('data/log.txt','a');
fprintf(file,'1D and 2D, both full.\n');
fprintf(file,'Elapsed time %02d:%02d:%02d.%03d\n',hr,min,sec,secdec);
fprintf(file,'data is a matrix with 1D in 1st col and 2D in 2nd col.\n\n\n');
fclose(file);

save('data/D12_full.mat','L','C','data')

%for j=1:n;fprintf('c = %3.1f \t\tDD = %1.8e\n',C(j),data(j));end


%% plots and asymptotics   (needs more love)
% Here we do see some numerical error in the eigenvalue calculation since
% they are not following the asymptotic formula given by Slepian IV (1964).
%pause()
%To be run later
%{
%C=linspace(0.1,16,100); 
%data=load('lambda_c0.1-16.tsv','-ascii');
plot(C,1-data)
hold on

%I=(20:50).';
%A=[C(I)', ones(size(I))]\log(1-data(I));
%plot(C, exp(A(2)+A(1)*C))

plot(C,pi*8*C.*exp(-2*C))%asymptotic formula eqn. (93) of Slepian IV (1964) 

set(gca,'xlim',[0,10],'yscale','log')
%}


%% D=2, not the full disc, c-loop 
clc;clear all

L=2000;
n=21;
%C=[.1 .5 1 1.5 2 3 5 10];
C=linspace(0,10,n);
R=1;

Q=[.1, .3, .5, .7, .9];
data=zeros(n,length(Q)); %init

file=fopen('data/log.txt','a');

for b=1:length(Q) 
clear A Gamma x

q=Q(b)
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
time=toc;
hr  = floor(time/3600); time=time-hr*3600;
min = floor(time/60); time=time-min*60;
sec=floor(time); 
secdec=floor(1000*(time-sec));
fprintf(file,'2D circle ring, q=%1.1f, calculating kernel matrix.\n', q);
fprintf(file,'Elapsed time %02d:%02d:%02d.%03d\n',hr,min,sec,secdec);

tic

parfor j=1:n
    Y=repmat(x{j}',L,1);
    data(j,b)=eigs(A{j}.*Y,1)*R^2*Gamma{j}*(1-q)/L;
end
time=toc;
hr  = floor(time/3600); time=time-hr*3600;
min = floor(time/60); time=time-min*60;
sec=floor(time); 
secdec=floor(1000*(time-sec));
fprintf(file,'2D eigenvalues, circle ring.\n');
fprintf(file,'Elapsed time %02d:%02d:%02d.%03d\n\n',hr,min,sec,secdec);



end

fprintf(file,'data is a matrix with each col corresponding to a q value in Q.\n\n\n');
fclose(file);
save('data/D2_ring.mat','L','C','Q','data')















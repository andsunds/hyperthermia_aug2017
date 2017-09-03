%% ring_T-H_L1000
clc;clf;clear all


load('data/D2_ring-bessel-AB-0_L-1000_2017-09-02.mat');
plot(C/2/pi,data,'k-')
hold on

%load('data/D2_ring-0_L-2000_2017-08-19.mat');
load('data/D2_ring-bessel-AB-0_L-1000_2017-09-01.mat');

plot(C/2/pi,data)

xtext='$R_\mathrm{T}/\lambda$';
xlabel(xtext, 'interpreter','latex')
ytext='$||E||_\mathrm{T}^2/||E||_\mathrm{H}^2$';
ylabel(ytext, 'interpreter','latex')
set(gca, 'fontsize',18, 'xlim',[0,2])

lcell=cellstr(num2str([0,Q]','$q=%0.1f$'));
l=legend(lcell);
set(l, 'location','northwest','interpreter','latex')

%% ring-both_L1000
clc;clf;clear all

load('data/D2_ring-bessel-0_L-1000_2017-09-02.mat')

plot(C/2/pi,data(:,1),'-k'), hold on
plot(C/2/pi,data(:,2:end))


xtext='$R_\mathrm{T}/\lambda_{\mathrm{min}}$';
xlabel(xtext, 'interpreter','latex')
ytext='Energy fraction';%'$||E||_\mathrm{T}^2/||E||_{\infty}^2$';
ylabel(ytext, 'interpreter','latex')
set(gca, 'fontsize',18, 'xlim',[0,2])

lcell=cellstr(num2str(Q','$q=%0.2f$'));
l=legend(lcell);
%set(l, 'location','NorthWest','interpreter','latex')
set(l, 'position',[0.195,.74,0,0],'interpreter','latex')

load('data/D2_ring-bessel-AB-0_L-1000_2017-09-02.mat');
plot(C/2/pi,data,'k--')
hold on

%load('data/D2_ring-0_L-2000_2017-08-19.mat');
load('data/D2_ring-bessel-AB-0_L-1000_2017-09-01.mat');
plot(0,0,0,0)
plot(C/2/pi,data,'--')





clear all; close all; clc;

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

addpath('TO_dyn')
addpath('FO_dyn')

load("2011TOs_cp.mat")
load("2011RKs_EO_cp.mat")

FOTO_prob=prob([prob.isTO]==1);
%%
% TO costate solutions
YY_TO=[TO_prob.y0];
LL_TO=YY_TO(8:14,:);

% FOTO costate solutions
YY_FOTO=[FOTO_prob.y0];
LL_FOTO=YY_FOTO(8:14,:);

% gamma_Limit
gL=zeros(size(TO_prob));
for i=1:length(gL)
    gL(i)=-1/max(TO_prob(i).S);
end

% gamma
gM=LL_FOTO./LL_TO;
g=mean(gM,1);
g_std=std(gM,1,1);

figure
subplot(3,1,1)
plot(et2MJD2000([TO_prob.t0]),gL,'--r')
hold on
plot(et2MJD2000([TO_prob.t0]),g,'-ob')
hold off
legend('$$\gamma_L$$','$$\gamma$$','interpreter','latex','location','best','fontsize',15)
grid on
grid minor
axis tight
subplot(3,1,2)
semilogy(et2MJD2000([TO_prob.t0]),g_std)
legend('$$\sigma_{\gamma}$$','interpreter','latex','location','best','fontsize',15)
grid on
grid minor
axis tight
subplot(3,1,3)
plot(et2MJD2000([TO_prob.t0]),g./gL,'-ob')
legend('$$\frac{\gamma}{\gamma_L}$$','interpreter','latex','location','best','fontsize',15)
grid on
grid minor
axis tight

ì%%
% FOTO constraint violation
df_FOTO=sqrt(sum([FOTO_prob.gamma].^2,1));

% multiples of gL
k=[1*(1+[-1e-1 -1e-5 -1e-10 -eps 0 eps 1e-10 1e-5 1e-1]) 2];

df_FOTO_test=zeros(length(k),length(gL));

for i=1:length(k)
    for j=1:length(gL)
        gamma_test=FO_ZFP(k(i)*gL(j)*LL_TO(:,j),FOTO_prob(j));
        df_FOTO_test(i,j)=norm(gamma_test);
    end
end
%%
lw=0.5*ones(size(k))+1.5*(k==1)+1.5*(k==2);
df_FOTO_test_m=mean(df_FOTO_test,2);
df_FOTO_test_std=std(df_FOTO_test,1,2);

figure
subplot(1,6,1)
semilogy(k,df_FOTO_test_m,'-o')
grid on
grid minor
ylabel('$$\overline{\|\Gamma_{k\gamma_L}\|}$$','interpreter','latex')
xlabel('$$k$$','interpreter','latex')

subplot(1,6,2)
semilogy(k,df_FOTO_test_std,'-o')
grid on
grid minor
ylabel('$$\sigma_{\|\Gamma_{k\gamma_L}\|}$$','interpreter','latex')
xlabel('$$k$$','interpreter','latex')

subplot(1,6,3:6)
for i=1:length(k)
    semilogy(et2MJD2000([TO_prob.t0]),df_FOTO_test(i,:),'-','LineWidth',lw(i))
    hold on
end
hold off
grid on
grid minor
axis tight
legendLabels=arrayfun(@(x) sprintf('k=1+%.6e', x-1),k,'UniformOutput', false);
legend(legendLabels,'Location','best','interpreter','latex','fontsize',10);
title('$${\|\Gamma_{k\gamma_L}\|}$$','interpreter','latex')
xlabel('$$t_0 [MJD2000]$$','interpreter','latex')

%%
idf1=find([df_FOTO_test(k==1,:)]>=1e-9);

prob_test=FOTO_prob(idf1(1));
prob_test.y0(8:14)=gL(idf1(1))*TO_prob(idf1(1)).y0(8:14);

DispRes(prob_test);
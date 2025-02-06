clear all; clc;

load("VP_res\VP_2012FO.mat")
load("VP_res\VP_2012TO.mat")

epsilon=0;

%--------------------------------------------------------------------------

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

% Non-skipped implicit TO solutions
ITO_id=find(strcmp({EO_prob.sts},'TO'));
ITO_prob=EO_prob(ITO_id);

% Corresponding TO solutions
TOr_prob=TO_prob(find(ismember([TO_prob.t0],[ITO_prob.t0])));

% EO sim
fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',50);

if epsilon==1
    for i=1:length(ITO_prob)
        ITO_prob(i).epsilon=epsilon;
        [llEO]=fsolve(@(ll) FO_ZFP(ll,ITO_prob(i)),ITO_prob(i).y0(8:14),fsopt);
        [~,~,ITO_prob(i)]=FO_ZFP(llEO,ITO_prob(i));
        ITO_prob(i)=DispRes(ITO_prob(i),0);
    end
end

gL=zeros(size(TO_prob));
for i=1:length(TO_prob)
    gL(i)=-1/max(TO_prob(i).S);
end

g=zeros(size(TOr_prob));
sig=zeros(size(TOr_prob));
for i=1:length(TOr_prob)
    g(i)=mean(ITO_prob(i).y0(8:14)./TOr_prob(i).y0(8:14));
    sig(i)=std(ITO_prob(i).y0(8:14)./TOr_prob(i).y0(8:14));
end

% Set default font size for axes labels, tick labels and legend
fs=15;
set(0,'DefaultAxesFontSize',fs);
set(0,'DefaultLegendFontSize',fs);

% Set the default linewidth
set(0,'DefaultLineLineWidth',2);

% Set the default text and legend interpreter
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

f=figure;
f.Position = [488,242,1200,400];
subplot(2,1,1)
semilogy(et2MJD2000([TO_prob.t0]),gL,'r')
hold on
semilogy(et2MJD2000([TO_prob.t0]),2*gL,'k--')
semilogy(et2MJD2000([TOr_prob.t0]),g,'bo-')
hold off
grid on
grid minor
xlim([min(et2MJD2000([TO_prob.t0])),max(et2MJD2000([TO_prob.t0]))])
ylabel('$\gamma$')
legend('$\gamma_L$','$2\gamma_L$','$\overline{\gamma}$','Location','northeast')

subplot(2,1,2)
semilogy(et2MJD2000([TOr_prob.t0]),sig,'bo-')
grid on
grid minor
xlim([min(et2MJD2000([TO_prob.t0])),max(et2MJD2000([TO_prob.t0]))])
xlabel('$t_0\,\left[MJD2000\right]$')
ylabel('$\sigma_\gamma$')

idM=find(sig==min(sig));
f=figure;
f.Position = [488,242,1200,400];
plot(ITO_prob(idM).tt,ITO_prob(idM).zz(:,8:14),'r')
hold on
plot(TOr_prob(idM).tt,g(idM).*TOr_prob(idM).zz(:,8:14),'b--','LineWidth',0.5)
hold off
grid on
grid minor
xlim([min(TOr_prob(idM).tt), max(TOr_prob(idM).tt)])
xlabel('$t\,\left[days\right]$')
ylabel('$\mathbf{\lambda}_{FO},\,\overline{\gamma}\mathbf{\lambda}_{TO}$')
title(sprintf('$\\sigma_\\gamma=%.4e$', min(sig)), 'Interpreter', 'latex')
legend('$\mathbf{\lambda}_{FO}\left(t\right)$','','','','','','','$\overline{\gamma}\mathbf{\lambda}_{TO}\left(t\right)$','interpreter','latex')

idm=find(sig==min(sig));
fprintf('Min sigma=%.4e\n',min(sig));
[~,~,ITO_prob(idm)]=FO_ZFP(g(idm).*TOr_prob(idm).y0(8:14),ITO_prob(idm));
ITO_prob(idm)=DispRes(ITO_prob(idm),0);
DispRes(ITO_prob(idm));

% User plot setting removal
set(0,'DefaultAxesFontSize','remove');
set(0,'DefaultLegendFontSize','remove');
set(0,'DefaultLineLineWidth','remove');
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultLegendInterpreter','remove');
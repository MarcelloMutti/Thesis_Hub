clear all; clc;

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

% 2000
load("CP_res\CP_2000FO.mat")
load("CP_res\CP_2000TO.mat")

% Non-skipped implicit TO solutions
ITO_id=find(strcmp({EO_prob.sts},'TO'));
ITO_prob=EO_prob(ITO_id);

% Corresponding TO solutions
TOr_prob=TO_prob(find(ismember([TO_prob.t0],[ITO_prob.t0])));

t01=et2MJD2000([TOr_prob.t0]);
DF1=zeros(3,length(ITO_prob));
for i=1:length(ITO_prob)
    g=mean(ITO_prob(i).y0(8:14)./TOr_prob(i).y0(8:14));
    df=FO_ZFP(g*TOr_prob(i).y0(8:14),ITO_prob(i));
    DF1(1,i)=norm(df(1:3))*LU;
    DF1(2,i)=norm(df(4:6))*LU/TU;
    DF1(3,i)=abs(df(7));
end

%% 2010
load("CP_res\CP_2010FO.mat")
load("CP_res\CP_2010TO.mat")

% Non-skipped implicit TO solutions
ITO_id=find(strcmp({EO_prob.sts},'TO'));
ITO_prob=EO_prob(ITO_id);

% Corresponding TO solutions
TOr_prob=TO_prob(find(ismember([TO_prob.t0],[ITO_prob.t0])));

t02=et2MJD2000([TOr_prob.t0]);
DF2=zeros(3,length(ITO_prob));
for i=1:length(ITO_prob)
    g=mean(ITO_prob(i).y0(8:14)./TOr_prob(i).y0(8:14));
    df=FO_ZFP(g*TOr_prob(i).y0(8:14),ITO_prob(i));
    DF2(1,i)=norm(df(1:3))*LU;
    DF2(2,i)=norm(df(4:6))*LU/TU;
    DF2(3,i)=abs(df(7));
end

%% 2011
load("CP_res\CP_2011FO.mat")
load("CP_res\CP_2011TO.mat")

% Non-skipped implicit TO solutions
ITO_id=find(strcmp({EO_prob.sts},'TO'));
ITO_prob=EO_prob(ITO_id);

% Corresponding TO solutions
TOr_prob=TO_prob(find(ismember([TO_prob.t0],[ITO_prob.t0])));

t03=et2MJD2000([TOr_prob.t0]);
DF3=zeros(3,length(ITO_prob));
for i=1:length(ITO_prob)
    g=mean(ITO_prob(i).y0(8:14)./TOr_prob(i).y0(8:14));
    df=FO_ZFP(g*TOr_prob(i).y0(8:14),ITO_prob(i));
    DF3(1,i)=norm(df(1:3))*LU;
    DF3(2,i)=norm(df(4:6))*LU/TU;
    DF3(3,i)=abs(df(7));
end

%% 2012
load("CP_res\CP_2012FO.mat")
load("CP_res\CP_2012TO.mat")

% Non-skipped implicit TO solutions
ITO_id=find(strcmp({EO_prob.sts},'TO'));
ITO_prob=EO_prob(ITO_id);

% Corresponding TO solutions
TOr_prob=TO_prob(find(ismember([TO_prob.t0],[ITO_prob.t0])));

t04=et2MJD2000([TOr_prob.t0]);
DF4=zeros(3,length(ITO_prob));
for i=1:length(ITO_prob)
    g=mean(ITO_prob(i).y0(8:14)./TOr_prob(i).y0(8:14));
    df=FO_ZFP(g*TOr_prob(i).y0(8:14),ITO_prob(i));
    DF4(1,i)=norm(df(1:3))*LU;
    DF4(2,i)=norm(df(4:6))*LU/TU;
    DF4(3,i)=abs(df(7));
end

%% 2014
load("CP_res\CP_2014FO.mat")
load("CP_res\CP_2014TO.mat")

% Non-skipped implicit TO solutions
ITO_id=find(strcmp({EO_prob.sts},'TO'));
ITO_prob=EO_prob(ITO_id);

% Corresponding TO solutions
TOr_prob=TO_prob(find(ismember([TO_prob.t0],[ITO_prob.t0])));

t05=et2MJD2000([TOr_prob.t0]);
DF5=zeros(3,length(ITO_prob));
for i=1:length(ITO_prob)
    g=mean(ITO_prob(i).y0(8:14)./TOr_prob(i).y0(8:14));
    df=FO_ZFP(g*TOr_prob(i).y0(8:14),ITO_prob(i));
    DF5(1,i)=norm(df(1:3))*LU;
    DF5(2,i)=norm(df(4:6))*LU/TU;
    DF5(3,i)=abs(df(7));
end

%%

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
f.Position=[0,0, 1200, 800];
subplot(3,1,1)
semilogy(t01,DF1(1,:))
hold on
semilogy(t02,DF2(1,:))
semilogy(t03,DF3(1,:))
semilogy(t04,DF4(1,:))
semilogy(t05,DF5(1,:))
hold off
grid on
grid minor
xlim([min([t01 t02 t03 t04 t05]) max([t01 t02 t03 t04 t05])])
% yticks(linspace(min([DF1(1,:) DF2(1,:) DF3(1,:) DF4(1,:) DF5(1,:)]), max([DF1(1,:) DF2(1,:) DF3(1,:) DF4(1,:) DF5(1,:)]), 3))
ylabel('$\|\Delta\mathbf{r}\left(t_f\right)\|,\,\left[km\right]$')

subplot(3,1,2)
semilogy(t01,DF1(2,:))
hold on
semilogy(t02,DF2(2,:))
semilogy(t03,DF3(2,:))
semilogy(t04,DF4(2,:))
semilogy(t05,DF5(2,:))
hold off
grid on
grid minor
xlim([min([t01 t02 t03 t04 t05]) max([t01 t02 t03 t04 t05])])
ylabel('$\|\Delta\mathbf{v}\left(t_f\right)\|,\,\left[km/s\right]$')

subplot(3,1,3)
semilogy(t01,DF1(3,:))
hold on
semilogy(t02,DF2(3,:))
semilogy(t03,DF3(3,:))
semilogy(t04,DF4(3,:))
semilogy(t05,DF5(3,:))
hold off
grid on
grid minor
xlim([min([t01 t02 t03 t04 t05]) max([t01 t02 t03 t04 t05])])
ylabel('$\left|\Delta {\lambda}_m\left(t_f\right)\right|$')
xlabel('$t_0,\,\left[MJD2000\right]$')

% User plot setting removal
set(0,'DefaultAxesFontSize','remove');
set(0,'DefaultLegendFontSize','remove');
set(0,'DefaultLineLineWidth','remove');
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultLegendInterpreter','remove');

fprintf('dr range %.4e %.4e\n',min([DF1(1,:) DF2(1,:) DF3(1,:) DF4(1,:) DF5(1,:)]),max([DF1(1,:) DF2(1,:) DF3(1,:) DF4(1,:) DF5(1,:)]))
fprintf('dv range %.4e %.4e\n',min([DF1(2,:) DF2(2,:) DF3(2,:) DF4(2,:) DF5(2,:)]),max([DF1(2,:) DF2(2,:) DF3(2,:) DF4(2,:) DF5(2,:)]))
fprintf('dlm range %.4e %.4e\n',min([DF1(3,:) DF2(3,:) DF3(3,:) DF4(3,:) DF5(3,:)]),max([DF1(3,:) DF2(3,:) DF3(3,:) DF4(3,:) DF5(3,:)]))
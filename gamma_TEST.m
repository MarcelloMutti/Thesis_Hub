clear all; clc;

% test last solution, introuce isapprox instead of isequal

load("CP_res\CP_2011FO.mat")
load("CP_res\CP_2011TO.mat")

epsilon=1;

%--------------------------------------------------------------------------

addpath("TO_dyn")
addpath("FO_dyn")

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',50);

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

ITO_id=find(strcmp({EO_prob.sts},'TO'));
ITO_prob=EO_prob(ITO_id);

% Corresponding TO solutions
TOr_prob=TO_prob(find(ismember([TO_prob.t0],[ITO_prob.t0])));

kk=zeros(size(TOr_prob));
mS=zeros(size(TOr_prob));

if epsilon==0
    wb1=waitbar(0,sprintf('FO(TO) testing of %.0f/%.0f',1,length(TOr_prob)));
else
    wb1=waitbar(0,sprintf('EO(TO) testing of %.0f/%.0f',1,length(TOr_prob)));
end

for i=1:length(TOr_prob)

    if epsilon==0
        wb1=waitbar(i/length(TOr_prob),wb1,sprintf('FO(TO) testing of %.0f/%.0f',i,length(TOr_prob)));
    else
        wb1=waitbar(i/length(TOr_prob),wb1,sprintf('EO(TO) testing of %.0f/%.0f',i,length(TOr_prob)));
    end

    gL=-1/max(TOr_prob(i).S);

    ll_TO=TOr_prob(i).y0(8:14);

    ll_g=(1+epsilon)*gL*ll_TO;

    ITO_prob(i).epsilon=epsilon;

    % [~,df,ex_flag]=fsolve(@(ll) FO_ZFP(ll,ITO_prob(i)),ll_g,fsopt);
    df_i=FO_ZFP(ll_g,ITO_prob(i));

    c=norm(df_i(1:3)*LU)<10 && norm(df_i(4:6)*LU/TU)<1e-3 && df_i(7)<1e-12;

    % if ex_flag<=0 || ~isequal(df,df_i)
    if ~c

        k=(1+epsilon)*2;
        kt=(1+epsilon);

        ex_flag=0;

        k1=(kt+(k-kt)/2);
        k2=k;

        c=0;

        while ex_flag<=0

            ll_g1=k1*gL*ll_TO;
            ll_g2=k2*gL*ll_TO;

            % [~,df_1,ex_flag1]=fsolve(@(ll) FO_ZFP(ll,ITO_prob(i)),ll_g1,fsopt);
            df_i1=FO_ZFP(ll_g1,ITO_prob(i));

            % [~,df_2,ex_flag2]=fsolve(@(ll) FO_ZFP(ll,ITO_prob(i)),ll_g2,fsopt);
            df_i2=FO_ZFP(ll_g2,ITO_prob(i));

            % c1=(ex_flag1>0) && isequal(df_1,df_i1);
            % c2=(ex_flag2>0) && isequal(df_2,df_i2);

            c1=norm(df_i1(1:3)*LU)<10 && norm(df_i1(4:6)*LU/TU)<1e-3 && df_i1(7)<1e-12;
            c2=norm(df_i2(1:3)*LU)<10 && norm(df_i2(4:6)*LU/TU)<1e-3 && df_i2(7)<1e-12;

            if ~c1 && ~c2 && ~c && k2<=10*kt% no convergence yet
                k2=kt+(k2-kt)*2;
                k1=kt+(k2-kt)/2;
            elseif k2>10*kt % no convergence at all
                kk(i)=NaN;
                ex_flag=1;
            elseif (~c1 && ~c2 && c) || k2-kt<10*eps % stop convergence
                kk(i)=k_old;
                ex_flag=1;
            elseif ~c1 && c2 % mid 1-2
                c=1;
                k_old=k2;
                k2=k1+(k2-k1)/2;
                k1=kt+(k2-kt)/2;
                ex_flag=0;
            else
                c=1;
                k_old=k1;
                k2=k1;
                k1=kt+(k2-kt)/2;
                ex_flag=0;
            end

        end

    else

        kk(i)=1+epsilon;

    end

    if ~isnan(kk(i))
        gg=kk(i)*gL*ll_TO;
        [~,~,ITO_prob(i)]=FO_ZFP(gg,ITO_prob(i));
        ITO_prob(i)=DispRes(ITO_prob(i),0);
        mS(i)=max(ITO_prob(i).S);
    else
        mS(i)=NaN;
    end

end

close(wb1);

% Set default font size for axes labels, tick labels and legend
fs=15;
set(0,'DefaultAxesFontSize',fs);
set(0,'DefaultLegendFontSize',fs);

% Set the default linewidth
set(0,'DefaultLineLineWidth',1);

% Set the default text and legend interpreter
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

f=figure;
f.Position = [488,242,1200,400];
subplot(2,1,1)
semilogy(et2MJD2000([TOr_prob(~isnan(kk)).t0]),kk(~isnan(kk)),'bo-')
hold on
semilogy(et2MJD2000([TOr_prob(~isnan(kk)).t0]),(1+epsilon)*ones(size(kk(~isnan(kk)))),'k--')
hold off
grid on
grid minor
xlim([min(et2MJD2000([TOr_prob(~isnan(kk)).t0])) max(et2MJD2000([TOr_prob(~isnan(kk)).t0]))])
if epsilon==0
    ylabel('$\gamma_{FO}/\gamma_L$')
else
    ylabel('$\gamma_{EO}/\gamma_L$')
end

subplot(2,1,2)
semilogy(et2MJD2000([TOr_prob(~isnan(kk)).t0]),abs(mS(~isnan(kk))+epsilon),'bo-')
grid on
grid minor
xlim([min(et2MJD2000([TOr_prob(~isnan(kk)).t0])) max(et2MJD2000([TOr_prob(~isnan(kk)).t0]))])
if epsilon==0
    ylabel('$\left|\max\left(S_{LEFO,0}\right)\right|$')
else
    ylabel('$\left|\max\left(S_{LEFO,1}\right)+1\right|$')
end
xlabel('$t_0\,\left[MJD2000\right]$')

% User plot setting removal
set(0,'DefaultAxesFontSize','remove');
set(0,'DefaultLegendFontSize','remove');
set(0,'DefaultLineLineWidth','remove');
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultLegendInterpreter','remove');
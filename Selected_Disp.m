clear all; clc;

load("CP_res\2014FO.mat")

% 8650 or 8900
t0=8900;

% 700 or []
tf=700;

%--------------------------------------------------------------------------

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

addpath("TO_dyn")
addpath("FO_dyn")

t0=MJD20002et(t0);
tf_ad=tf*86400/TU;

fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-8,'FunctionTolerance',1e-8,'MaxIterations',2e2);

if isempty(tf)
    prob=TO_prob;
    id=find(abs([prob.t0]-t0)==min(abs(([prob.t0]-t0))),1,'first');

    sol=prob(id);
    sol.t0=t0;

    if sol.t0>prob(id).t0
        lltf_g=[prob(id).y0(8:14)+(prob(id+1).y0(8:14)-prob(id).y0(8:14))./2; prob(id).tf_ad+(prob(id+1).tf_ad-prob(id).tf_ad)/2];
    else
        lltf_g=[prob(id-1).y0(8:14)+(prob(id).y0(8:14)-prob(id-1).y0(8:14))./2; prob(id-1).tf_ad+(prob(id).tf_ad-prob(id-1).tf_ad)/2];
    end
    
    [lltf_TO,df]=fsolve(@(llt) TO_ZFP(llt,sol),lltf_g,fsopt);
    [~,~,sol]=TO_ZFP(lltf_TO,sol);
    sol=DispRes(sol);

else

    prob=EO_prob;
    idtf=find(abs([prob.tf_ad]-tf_ad)==min(abs(([prob.tf_ad]-tf_ad))));
    idt0=find(abs([prob(idtf).t0]-t0)==min(abs(([prob(idtf).t0]-t0))));
    id=idtf(idt0);

    sol=prob(id);
    sol.t0=t0;
    sol.tf_ad=tf_ad;
    sol.tf=t0+tf_ad*TU;
    sol.epsilon=1;
    ll_g=prob(id).y0(8:14);

    % if sol.t0>prob(id).t0
    %     ll_g=prob(id).y0(8:14)+(prob(idtf(idt0+1)).y0(8:14)-prob(id).y0(8:14))./2;
    % else
    %     ll_g=prob(idtf(idt0-1)).y0(8:14)+(prob(id).y0(8:14)-prob(idtf(idt0-1)).y0(8:14))./2;
    % end

    [ll_TO,df]=fsolve(@(llt) FO_ZFP(llt,sol),ll_g,fsopt);
    [~,~,sol]=FO_ZFP(ll_TO,sol);
    sol=DispRes(sol,0);

    sol=E2F_CONT(sol,1,1);

    DispRes(sol);
end




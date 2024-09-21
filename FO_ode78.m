function [tt,zz]=FO_ode78(prob,tspan,z0)

    epsilon=prob.epsilon;

    r0=norm(z0(1:3));
    [~,~,Sp]=MARGO_param(r0);

    if Sp<prob.Plim(2)
        Ptype='med';
    elseif Sp>=prob.Plim(2)
        Ptype='max';
    end

    if SwFun(tspan(1),z0,prob.isFO)+epsilon<0
        utype='on';
    elseif SwFun(tspan(1),z0,prob.isFO)-epsilon>0 || (SwFun(tspan(1),z0,prob.isFO)-epsilon>=0 && epsilon==0)
        utype='off';
    elseif abs(SwFun(tspan(1),z0,prob.isFO))-epsilon<=0 && epsilon~=0
        utype='med';
    end

    tt=[];
    zz=[];
    complete=0;
    ie_old=[];

    while ~complete

        options=odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@(t,y) crossings(t,y,prob,epsilon,ie_old));

        [tti,zzi,te,ze,ie]=ode78(@(t,z) FO_2BP_SEP(t,z,Ptype,utype,epsilon),tspan,z0,options);

        s_id=find(ie~=ie_old,1,'last');
        if isempty(s_id)
            s_id=1;
        end
    
        if ~isempty(ie)

            tt=[tt; tti(1:end)];   % to try correcting for duplicates
            zz=[zz; zzi(1:end,:)];

            ti=te(s_id);   % possibly picks up 0 at t=0
            ze=ze(s_id,:);

            z0=ze.';
            tspan(1)=ti;

            switch ie(s_id)

                case 1*~strcmp(Ptype,'max')  % Pmed -> Pmax

                    re=ze(1:3).';
                    ve=ze(4:6).';

                    dze_m=FO_2BP_SEP(te,ze.','med',utype,epsilon);
                    dze_p=FO_2BP_SEP(te,ze.','max',utype,epsilon);

                    Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*[re.', zeros(1,11)]./(re.'*ve);

                    Phie_m=reshape(ze(15:210),[14,14]);
                    Phie_p=Psi*Phie_m;
                    z0(15:210)=reshape(Phie_p,[1,14*14]);

                    Ptype='max';

                case 2*~strcmp(Ptype,'med')  % Pmax -> Pmed

                    re=ze(1:3).';
                    ve=ze(4:6).';

                    dze_m=FO_2BP_SEP(te,ze.','max',utype,epsilon);
                    dze_p=FO_2BP_SEP(te,ze.','med',utype,epsilon);

                    Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*[re.', zeros(1,11)]./(re.'*ve); 

                    Phie_m=reshape(ze(15:210),[14,14]);
                    Phie_p=Psi*Phie_m;
                    z0(15:210)=reshape(Phie_p,[1,14*14]);

                    Ptype='med';

                case 3*~strcmp(utype,'off')*~strcmp(utype,'med')  % u_on -> u_med eps~0 (u_on -> u_off eps=0)

                    if epsilon==0   % u_on -> u_off

                        re=ze(1:3).';
                        ve=ze(4:6).';
                        me=ze(7);
                        lre=ze(8:10).';
                        lve=ze(11:13).';
    
                        dze_m=FO_2BP_SEP(te,ze.',Ptype,'on',epsilon);
                        dze_p=FO_2BP_SEP(te,ze.',Ptype,'off',epsilon);

                        [Tc,Tcp]=MARGO_param(norm(re));

                        c=Tc(2);
                        cp=Tcp(2);

                        DyS=[-norm(lve)*cp/me*re.'/norm(re), zeros(1,3), norm(lve)*c/me^2, zeros(1,3), -c/me*lve.'/norm(lve), -1];
                        DtS=norm(lve)*c/me*lre.'*lve-norm(lve)/me*cp*re.'*ve/norm(re);
                        
    
                        Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*DyS/DtS;
    
                        Phie_m=reshape(ze(15:210),[14,14]);
                        Phie_p=Psi*Phie_m;
                        z0(15:210)=reshape(Phie_p,[1,14*14]);

                        utype='off';

                    else    % u_on -> u_med

                        utype='med';

                    end

                case 4*(~strcmp(utype,'on')*~strcmp(utype,'off')*(epsilon>0)+~strcmp(utype,'on')*~strcmp(utype,'med')*(epsilon==0))  % u_med -> u_on eps~0 (u_off -> u_on eps=0)

                    if epsilon==0   % u_off -> u_on

                        re=ze(1:3).';
                        ve=ze(4:6).';
                        me=ze(7);
                        lre=ze(8:10).';
                        lve=ze(11:13).';
    
                        dze_m=FO_2BP_SEP(te,ze.',Ptype,'off',epsilon);
                        dze_p=FO_2BP_SEP(te,ze.',Ptype,'on',epsilon);

                        [Tc,Tcp]=MARGO_param(norm(re));

                        c=Tc(2);
                        cp=Tcp(2);

                        DyS=[-norm(lve)*cp/me*re.'/norm(re), zeros(1,3), norm(lve)*c/me^2, zeros(1,3), -c/me*lve.'/norm(lve), -1];
                        DtS=norm(lve)*c/me*lre.'*lve-norm(lve)/me*cp*re.'*ve/norm(re);
                        
    
                        Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*DyS/DtS;
    
                        Phie_m=reshape(ze(15:210),[14,14]);
                        Phie_p=Psi*Phie_m;
                        z0(15:210)=reshape(Phie_p,[1,14*14]);

                        utype='on';

                    else    % u_med -> u_on

                        utype='on';

                    end

                case 5*~strcmp(utype,'on')*~strcmp(utype,'off')  % u_med -> u_off

                    utype='off';

                case 6*~strcmp(utype,'on')*~strcmp(utype,'med')  % u_off -> u_med

                    utype='med';

                otherwise
            end

            ie_old=ie(s_id);

        else

            tt=[tt; tti];
            zz=[zz; zzi];
            complete=1;

        end

    end
    
end

function [value,isterminal,direction]=crossings(t,y,prob,epsilon,ie_old)

    rr=y(1:3);
    vv=y(4:6);
    m=y(7);
    llr=y(8:10);
    llv=y(11:13);

    r=norm(rr);
    [Tc,Tcp,Sp,dSpdr]=MARGO_param(r);

    dSp=dSpdr*dot(rr,vv)/norm(rr);

    f1=Sp-prob.Plim(2);

    c=Tc(2);
    cp=Tcp(2);

    dSe=c/(m*norm(llv))*dot(llr,llv)-norm(llv)/m*cp*dot(rr,vv)/norm(rr);

    f2=SwFun(t,y,prob.isFO);

    if epsilon>0

        value=[f1/abs(dSp);  % Power switching function                      
               f1/abs(dSp);
               (f2+epsilon)/abs(dSe);  % throttle switching function (for TO)
               (f2+epsilon)/abs(dSe);
               (f2-epsilon)/abs(dSe);
               (f2-epsilon)/abs(dSe)];
    
        isterminal=ones(6,1);
        
        direction=[+1;  % from Pmed to Pmax [ie=1]
                   -1;  % from Pmax to Pmed [ie=2]
                   +1;  % on to med         [ie=3]
                   -1;  % med to on         [ie=4]
                   +1;  % med to off        [ie=5]
                   -1]; % off to med        [ie=6]

    else

        value=[f1/dSp;  % Power switching function                      
               f1/dSp;
               (f2+epsilon)/dSe;  % throttle switching function (for TO)
               (f2+epsilon)/dSe];
    
        isterminal=ones(4,1);
        
        direction=[+1;  % from Pmed to Pmax [ie=1]
                   -1;  % from Pmax to Pmed [ie=2]
                   +1;  % on to off         [ie=3]
                   -1]; % off to on         [ie=4]

    end

    value(ie_old)=inf;

end
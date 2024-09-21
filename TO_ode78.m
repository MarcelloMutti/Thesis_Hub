function [tt,zz]=TO_ode78(prob,tspan,z0)

    r0=norm(z0(1:3));
    [~,~,Sp]=MARGO_param(r0);

    if Sp<prob.Plim(2)
        Ptype='med';
    elseif Sp>=prob.Plim(2)
        Ptype='max';
    end

    tt=[];
    zz=[];
    complete=0;
    ie_old=[];

    while ~complete

        options=odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@(t,y) crossings(t,y,prob,ie_old));

        [tti,zzi,te,ze,ie]=ode78(@(t,z) TO_2BP_SEP(t,z,Ptype),tspan,z0,options);

        s_id=find(ie~=ie_old,1,'last');
        if isempty(s_id)
            s_id=1;
        end
    
        if ~isempty(ie)

            tt=[tt; tti(1:end)];   % to try correcting for duplicates
            zz=[zz; zzi(1:end,:)];

            ti=te(s_id);
            ze=ze(s_id,:);

            z0=ze.';
            tspan(1)=ti;

            switch ie(s_id)

                case 1*~strcmp(Ptype,'max')  % Pmed -> Pmax

                    re=ze(1:3).';
                    ve=ze(4:6).';

                    dze_m=TO_2BP_SEP(te,ze.','med');
                    dze_p=TO_2BP_SEP(te,ze.','max');

                    Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*[re.', zeros(1,11)]./(re.'*ve);

                    Phie_m=reshape(ze(15:210),[14,14]);
                    Phie_p=Psi*Phie_m;
                    z0(15:210)=reshape(Phie_p,[1,14*14]);

                    Ptype='max';

                case 2*~strcmp(Ptype,'med')  % Pmax -> Pmed

                    re=ze(1:3).';
                    ve=ze(4:6).';

                    dze_m=TO_2BP_SEP(te,ze.','max');
                    dze_p=TO_2BP_SEP(te,ze.','med');

                    Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*[re.', zeros(1,11)]./(re.'*ve);

                    Phie_m=reshape(ze(15:210),[14,14]);
                    Phie_p=Psi*Phie_m;
                    z0(15:210)=reshape(Phie_p,[1,14*14]);

                    Ptype='med';

%                 case 3  % u_on -> u_off
%                 case 4  % u_off -> u_on
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

function [value,isterminal,direction]=crossings(t,y,prob,ie_old)

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

    value=[f1/abs(dSp);  % Power switching function                      
           f1/abs(dSp);
           f2/abs(dSe);  % throttle switching function (for FO)
           f2/abs(dSe)];

    isterminal=ones(4,1);
    
    direction=[+1;  % from Pmed to Pmax [ie=1]
               -1;  % from Pmax to Pmed [ie=2]
               +1;  % on to off         [ie=3]
               -1]; % off to on         [ie=4]

    value(ie_old)=inf;

end
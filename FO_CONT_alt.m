function [prob]=FO_CONT_alt(prob,TO_ref)
    
    % continuate along each column
    L=length(TO_ref);
    Lp=0;

    for id=1:L

        prob=EO_tfCONT_alt(prob,TO_ref,id);
        prob(Lp+1).isTO=1;
        % prob(Lp+1)=E2F_CONT(prob(Lp+1),id,L);
        Lp=length(prob);

    end

    % for i=1:length(prob)
    %     prob(i)=E2F_CONT(prob(i),i,length(prob));
    % end

end


function [prob]=EO_tfCONT_alt(prob,TO_ref,id)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    it=1;
    L=length(prob);
    if id==1
        L=0;
    end

    Dt_max=20; % [days]
    Dt_min=1;  % [days]
    ToF_max=800;

    iscomplete=0;
    % skip=0;
    nang=0;

    wb1=waitbar(0,sprintf('EO tf continuation [%.2f %%] of %.0f/%.0f',0,id,length(TO_ref)));

    while ~iscomplete

        prob(L+it)=prob(end);

        if it==1        %-analytical-rel-exploit---------------------------

            ex_flag=0;
            f=1;

            prob(L+it).tf_ad=TO_ref(id).tf_ad;
            prob(L+it).tf=TO_ref(id).tf;
            prob(L+it).t0=TO_ref(id).t0;
            % prob(L+it).epsilon=0;   % direct FO cont imp; improper to cont on EO from init FO sol;

            while ex_flag<=0

                gL=-1/max(TO_ref(id).S);
                % k=3/2+1/2*tanh((f/10-1/2)/0.1);
                k=(prob(L+it).epsilon+1)*1.01^(f-1);

                ll_g=k*gL*TO_ref(id).y0(8:14);

                [ll_FO,df,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);

                % df=FO_ZFP(ll_FO,prob(L+it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                

                % f=f+1;

            end

        elseif it==2

            prob(L+it).epsilon=1;

            ex_flag=0;
            f=1;

            while ex_flag<=0

                prob(L+it).tf_ad=min(ToF_max*86400/TU,prob(L+it-1).tf_ad+max(Dt_min*86400/TU,DT/(2^(f-1))));
                prob(L+it).tf=prob(L+it).tf_ad*TU+prob(L+it).t0;

                if f==1 && ~nang % attempt 0NPCM

                    ll_g=prob(L+it-1).y0(8:14);

                else    % ACT

                    ll_g=ACT(prob(L+it));

                end

                if anynan(ll_g)
                    nang=1;
                    ex_flag=0;
                else
                    [ll_FO,df,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);       
                end                

                % df=FO_ZFP(ll_FO,prob(L+it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                % if prob(L+it).tf_ad-TO_ref(id).tf_ad<Dt_max*86400/TU && f>10
                %     ex_flag=1;
                %     skip=1;
                % end
                

                f=f+1;

            end

        else    %-1-3NPCM--------------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0

                prob(L+it).tf_ad=min(ToF_max*86400/TU,prob(L+it-1).tf_ad+max(Dt_min*86400/TU,DT/(2^(f-1))));
                prob(L+it).tf=prob(L+it).tf_ad*TU+prob(L+it).t0;

                if rem(f,4-(it<3)-(it<4))==1 && ~nang % attempt 1NPCM

                    ll_g=prob(L+it-1).y0(8:14)+(prob(L+it).tf_ad-prob(L+it-1).tf_ad)*(prob(L+it-1).y0(8:14)-prob(L+it-2).y0(8:14))/(prob(L+it-1).tf_ad-prob(L+it-2).tf_ad);            
                   
                elseif rem(f,4-(it<4))==2 && it>=3 && ~nang % attempt 2NPCM

                    yy=[prob(L+it-3:L+it-1).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(L+it-3:L+it-1).tf_ad],ll,prob(L+it).tf_ad);

                elseif rem(f,4)==3 && it>=4 && ~nang % attempt 3NPCM

                    yy=[prob(L+it-4:L+it-1).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(L+it-4:L+it-1).tf_ad],ll,prob(it).tf_ad);

                else    % ACT

                    ll_g=ACT(prob(L+it));

                end

                if anynan(ll_g)
                    nang=1;
                    ex_flag=0;
                else
                    [ll_FO,df,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);       
                end                

                % df=FO_ZFP(ll_FO,prob(L+it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                % if prob(L+it).tf_ad-TO_ref(id).tf_ad<Dt_min*86400/TU && f>10
                %     ex_flag=1;
                %     skip=1;
                % end
                

                f=f+1;

            end

        end


        [~,~,prob(L+it)]=FO_ZFP(ll_FO,prob(L+it));

        prob(L+it)=DispRes(prob(L+it),0);

%         prob(L+it+1)=prob(L+it);
        if it==1
            DT=Dt_min*86400/TU;
        else
            DT=max(min(1.25*(prob(L+it).tf_ad-prob(L+it-1).tf_ad),Dt_max*86400/TU),Dt_min*86400/TU);
        end

        nang=0;

        if prob(L+it).tf_ad*TU/86400==ToF_max

            iscomplete=1;

        end

        % if skip
        % 
        %     iscomplete=1;
        %     prob=prob(1:end-1);
        %     prob(end).sts='skp';
        % 
        % end

        
        wb1=waitbar((prob(L+1).tf_ad-prob(L+it).tf_ad)/(prob(L+1).tf_ad-ToF_max*86400/TU),wb1,sprintf('EO tf continuation [%.2f %%] of %.0f/%.0f',abs(prob(L+1).tf_ad-prob(L+it).tf_ad)/abs(prob(L+1).tf_ad-ToF_max*86400/TU)*100,id,length(TO_ref)));

        it=it+1;     
        
    end


    prob(L+1).sts='TO';                 % obsolete

    fprintf('\n')

    close(wb1);

end


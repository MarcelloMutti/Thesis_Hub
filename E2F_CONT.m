function prob = E2F_CONT(prob)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-8,'FunctionTolerance',1e-8,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    it=1;

    E_min=1e-6;
    DE=0.05;
%     E=prob.epsilon;

    iscomplete=0;

    while ~iscomplete

        prob(it+1)=prob(it);

        if it==1        %-0NPCM--------------------------------------------

            ex_flag=0;
            f=1;
            E=prob(it).epsilon;

            while ex_flag<=0

                prob(it+1).epsilon=max(E-DE/(2^(f-1)),E_min);

                if rem(f,2)==1 % attempt 0NPCM (alternates with act)

                    ll_g=prob(it).y0(8:14);        
                   
                else    % ACT

                    ll_g=ACT(prob(it+1));

                end

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it+1)),ll_g,fsopt);                

                df=FO_ZFP(ll_FO,prob(it+1));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end

                f=f+1;

            end


             
%             wb2=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),'Initiating continuation');

%         elseif it==2    %-0NPCM--------------------------------------------
%             
%             tic
% 
%             ex_flag=0;
%             f=1;
% 
%             while ex_flag<=0
% 
%                 Dt_atmp=Dt_iter/(2^(f-1));
%                 prob(it).t0=prob(it-1).t0+Dt_atmp;
% 
%                 if rem(f,2)==1 % attempt 0NPCM (alternates with act)
% 
%                     ll_g=[prob(it-1).y0(8:14); prob(it-1).tf_ad];        
%                    
%                 else    % ACT
% 
%                     ll_g=[ACT(prob(it)); prob(it-1).tf_ad];
% 
%                 end
% 
%                 [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it)),ll_g,fsopt);
% 
%                 
% 
%                 df=FO_ZFP(ll_FO,prob(it));
% 
%                 if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
%                     ex_flag=0;
%                 end
% 
%                 f=f+1;
% 
%             end
% 
%         elseif it==2
% 
%             ex_flag=0;
%             f=1;
%             E=prob(it).epsilon;
% 
%             while ex_flag<=0
% 
%                 prob(it).epsilon=E-DE/(2^(f-1));
% 
%                 if rem(f,2)==1 % attempt 1NPCM
% 
%                     ll_g=prob(it-1).y0(8:14)+(prob(it).tf_ad-prob(it-1).tf_ad)*(prob(it-1).y0(8:14)-prob(id).y0(8:14))/(prob(it-1).tf_ad-prob(id).tf_ad);            
%                    
%                 else    % ACT
% 
%                     ll_g=ACT(prob(it));
% 
%                 end
% 
%                 [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it)),ll_g,fsopt);
%                 
% 
%                 df=FO_ZFP(ll_FO,prob(it));
% 
%                 if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
%                     ex_flag=0;
%                 end
%                 if prob(it).tf_ad-TO_ref(id).tf_ad<Dt_max*86400/TU && f>10
%                     ex_flag=1;
%                     skip=1;
%                 end
%                 
% 
%                 f=f+1;
% 
%             end

        else    %-1-3NPCM--------------------------------------------------

            ex_flag=0;
            f=1;
            E=prob(it).epsilon;

            while ex_flag<=0

                prob(it+1).epsilon=max(E-DE/(2^(f-1)),E_min);

                if rem(f,4)==1 % attempt 1NPCM

                    ll_g=prob(it).y0(8:14)+(prob(it+1).epsilon-prob(it).epsilon)*(prob(it).y0(8:14)-prob(it-1).y0(8:14))/(prob(it).epsilon-prob(it-1).epsilon);            
                   
                elseif rem(f,4)==2 && it>=3 % attempt 2NPCM

                    yy=[prob(it-2:it).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(it-2:it).epsilon],ll,prob(it+1).epsilon);

                elseif rem(f,4)==3 && it>=4 % attempt 3NPCM

                    yy=[prob(it-3:it).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(it-3:it).epsilon],ll,prob(it+1).epsilon);

                else    % ACT

                    ll_g=ACT(prob(it+1));

                end

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it+1)),ll_g,fsopt);
                

                df=FO_ZFP(ll_FO,prob(it+1));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end                

                f=f+1;

            end

        end


        [~,~,prob(it+1)]=FO_ZFP(ll_FO,prob(it+1));

        prob(it+1)=DispRes(prob(it+1),0);

        DE=1.25*(prob(it).epsilon-prob(it+1).epsilon);

%         prob(it+1)=prob(it);



        if prob(it+1).epsilon==E_min

            iscomplete=1;
            
        end


        it=it+1;

        


        
    end

    prob=prob(end);


    
    

end

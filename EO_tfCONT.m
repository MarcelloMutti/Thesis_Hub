function [prob]=EO_tfCONT(prob,TO_ref,id)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    it=1;
    L=length(prob);

    Dt_max=30; % [days]
    Dt_min=1;  % [days]

    iscomplete=0;
    skip=0;

    while ~iscomplete

        prob(L+it)=prob(id);

        if it==1        %-0NPCM--------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0
                
                prob(L+it).tf_ad=max(TO_ref(id).tf_ad,prob(id).tf_ad-(Dt_max/(2^(f-1)))*86400/TU);
                prob(L+it).tf=prob(L+it).tf_ad*TU+prob(L+it).t0;

                if rem(f,2)==1 % attempt 0NPCM (alternates with act)

                    ll_g=prob(id).y0(8:14);        
                   
                else    % ACT

                    ll_g=ACT(prob(L+it));

                end

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);                

                df=FO_ZFP(ll_FO,prob(L+it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                if prob(L+it).tf_ad-TO_ref(id).tf_ad<Dt_max*86400/TU && f>10
                    ex_flag=1;
                    skip=1;
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

        elseif it==2

            ex_flag=0;
            f=1;

            while ex_flag<=0

                prob(L+it).tf_ad=max(TO_ref(id).tf_ad,prob(L+it-1).tf_ad-DT/(2^(f-1)));
                prob(L+it).tf=prob(L+it).tf_ad*TU+prob(L+it).t0;

                if rem(f,2)==1 % attempt 1NPCM

                    ll_g=prob(L+it-1).y0(8:14)+(prob(L+it).tf_ad-prob(L+it-1).tf_ad)*(prob(L+it-1).y0(8:14)-prob(id).y0(8:14))/(prob(L+it-1).tf_ad-prob(id).tf_ad);            
                   
                else    % ACT

                    ll_g=ACT(prob(L+it));

                end

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);
                

                df=FO_ZFP(ll_FO,prob(L+it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                if prob(L+it).tf_ad-TO_ref(id).tf_ad<Dt_max*86400/TU && f>10
                    ex_flag=1;
                    skip=1;
                end
                

                f=f+1;

            end

        else    %-1-3NPCM--------------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0

                prob(L+it).tf_ad=max(TO_ref(id).tf_ad,prob(L+it-1).tf_ad-DT/(2^(f-1)));
                prob(L+it).tf=prob(L+it).tf_ad*TU+prob(L+it).t0;

                if rem(f,4)==1 % attempt 1NPCM

                    ll_g=prob(L+it-1).y0(8:14)+(prob(L+it).tf_ad-prob(L+it-1).tf_ad)*(prob(L+it-1).y0(8:14)-prob(L+it-2).y0(8:14))/(prob(L+it-1).tf_ad-prob(L+it-2).tf_ad);            
                   
                elseif rem(f,4)==2 && it>=3 % attempt 2NPCM

                    yy=[prob(L+it-3:L+it-1).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(L+it-3:L+it-1).tf_ad],ll,prob(L+it).tf_ad);

                elseif rem(f,4)==3 && it>=4 % attempt 3NPCM

                    yy=[prob(L+it-4:L+it-1).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(L+it-4:L+it-1).tf_ad],ll,prob(it).tf_ad);

                else    % ACT

                    ll_g=ACT(prob(L+it));

                end

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);
                

                df=FO_ZFP(ll_FO,prob(L+it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                if prob(L+it).tf_ad-TO_ref(id).tf_ad<Dt_max*86400/TU && f>10
                    ex_flag=1;
                    skip=1;
                end
                

                f=f+1;

            end

        end


        [~,~,prob(L+it)]=FO_ZFP(ll_FO,prob(L+it));

        prob(L+it)=DispRes(prob(L+it),0);

%         prob(L+it+1)=prob(L+it);
        if it==1
            DT=min(1.25*(prob(L+it).tf_ad-TO_ref(id).tf_ad),Dt_max*86400/TU);
        else
            DT=max(min(1.25*(prob(L+it-1).tf_ad-prob(L+it).tf_ad),Dt_max*86400/TU),Dt_min*86400/TU);
        end

        if prob(L+it).tf_ad==TO_ref(id).tf_ad

            iscomplete=1;
            prob(end).sts='TO';

        end

        if skip

            iscomplete=1;
            prob=prob(1:end-1);
            prob(end).sts='skp';

        end


        it=it+1;
            
%             if it~=1
%                 Dt_iter=Dt_atmp;
% 
%                 if (f-1)==1
%                     Dt_iter=min(1.1*Dt_iter,Dt_max*86400);
%                 end
%             end
% 
%             wb2=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),wb2,'TO continuation');
% 
%             it=it+1;
%             
%         end

        


        
    end

%     fprintf('\n')
% 
%     close(wb2);
% 
%     figure
%     plot(et2MJD2000([prob.t0]),[prob.tf_ad]*TU/86400,'linewidth',2)
%     grid on
%     grid minor
%     axis tight
%     ylim([100 1100])
%     title('tf')
% 
%     s=length(prob);
%     id=[1:10:s, s];
% 
%     figure
%     for i=id
%         plot3(prob(i).zz(:,1),prob(i).zz(:,2),prob(i).zz(:,3),'color',[.7 .7 .7])
%         view([55, 55])
%         hold on
%         plot3(prob(i).zz(1,1),prob(i).zz(1,2),prob(i).zz(1,3),'ob')
%         plot3(prob(i).zz(end,1),prob(i).zz(end,2),prob(i).zz(end,3),'kx')
%         plot3(0,0,0,'+k')
%     end
%     grid on
%     grid minor
%     xlabel('$x [AU]$')
%     ylabel('$y [AU]$')
%     zlabel('$z [AU]$')
% 
%     toc


    
    

end
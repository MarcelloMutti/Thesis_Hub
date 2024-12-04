function [prob]=FO_CONT(prob,TO_ref)
        
    % create row of max ToF
    prob=EO_t0CONT(prob,TO_ref);
    
    % continuate along each column
    L=length(prob);

    for id=1:L

        prob=EO_tfCONT(prob,TO_ref,id);
        prob(end).isTO=1;
        % prob(end)=E2F_CONT(prob(end),id,L);

    end

    for i=1:length(prob)
        prob(i)=E2F_CONT(prob(i),i,length(prob));
    end

end

function [prob] = EO_t0CONT(prob,TO_ref)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',2e2);
    % fsopt2=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

%     it=1;
% 
    t_wo=prob(1).tw(1);
    t_wc=prob(1).tw(2);
% 
%     prob(it).t0=t_wo;
% 
%     Dt_max=1; % [days]
%     Dt_iter=86400;
%     N=25;

    L=length(TO_ref);

    ToF_m=800; % [d]

    for it=1:L

        prob(it).t0=TO_ref(it).t0;
        prob(it).tf=TO_ref(it).t0+ToF_m*86400;
        prob(it).tf_ad=ToF_m*86400/TU;
%         prob(it).epsilon=1;

        if it==1        %-1st-solution-           
                        
            ex_flag=0;
            f=1;

            wb1=waitbar(0,'Generating first EO solution');

            while ex_flag<=0

                llg=ACT(prob(it));

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it)),llg,fsopt);

                df=FO_ZFP(ll_FO,prob(it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end

                f=f+1;

            end

            wb1=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),wb1,'Initiating EO t0 continuation');

        elseif it==2    %-0NPCM--------------------------------------------
            
            tic

            ex_flag=0;
            f=1;

            while ex_flag<=0

                if f==1 % attempt 0NPCM (alternates with act)

                    llg=prob(it-1).y0(8:14);        
                   
                else    % ACT

                    llg=ACT(prob(it));

                end

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it)),llg,fsopt);                

                df=FO_ZFP(ll_FO,prob(it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end

                f=f+1;

            end


        else    %-1-3NPCM--------------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0

                if f==1 % attempt 1NPCM

                    llg=prob(it-1).y0(8:14)+(prob(it).t0-prob(it-1).t0)*(prob(it-1).y0(8:14)-prob(it-2).y0(8:14))/(prob(it-1).t0-prob(it-2).t0);            
                   
                elseif f==2 && it>=4 % attempt 2NPCM

                    yy=[prob(it-3:it-1).y0];
                    ll=yy(8:14,:);
                    llg=makima([prob(it-3:it-1).t0],ll,prob(it).t0);

                elseif f==3 && it>=5 % attempt 3NPCM

                    yy=[prob(it-4:it-1).y0];
                    ll=yy(8:14,:);
                    llg=makima([prob(it-4:it-1).t0],ll,prob(it).t0);

                else    % ACT

                    llg=ACT(prob(it));

                end

                [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it)),llg,fsopt);
                

                df=FO_ZFP(ll_FO,prob(it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end

                f=f+1;

            end

        end

        if it~=L
            prob(it+1)=prob(it);
        end



        [~,~,prob(it)]=FO_ZFP(ll_FO,prob(it));

        prob(it)=DispRes(prob(it),0);

        wb1=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),wb1,sprintf('EO t0 continuation [%.2f %%]',(prob(it).t0-t_wo)/(t_wc-t_wo)*100));

%         if it==2
%             prob(1)=E2F_CONT(prob(1),1,1);
%             prob(2)=E2F_CONT(prob(2),1,1);
%         end
% 
% 
% 
% 
%         if it==1   % EO -> FO
% 
%             de=0.05;
%             
%             while prob(it).epsilon~=0
% 
%                 ex_flag=0;
%                 f=1;
% %                 de=0.05;
%                 ep=prob(it).epsilon;
%     
%                 while ex_flag<=0
% 
%                     prob(it).epsilon=max(ep-de/(2^(f-1)),0);
%     
%                     llg=prob(it).y0(8:14);
%         
%                     [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it)),llg,fsopt);
%                         
%         
%                     df=FO_ZFP(ll_FO,prob(it));
%         
%                     if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
%                         ex_flag=0;
%                     end
%         
%                     f=f+1;
%     
%                 end
% 
%                 [~,~,prob(it)]=FO_ZFP(ll_FO,prob(it));
% 
%                 prob(it)=DispRes(prob(it),0);
% 
%             end
%     
%             
%         end
% 
%         if prob(it).t0<t_wc
%             prob(it+1)=prob(it);
%             
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

    fprintf('\n')

    close(wb1);

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

end

function [prob]=EO_tfCONT(prob,TO_ref,id)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    it=1;
    L=length(prob);

    Dt_max=20; % [days]
    Dt_min=1;  % [days]

    iscomplete=0;
    skip=0;
    nang=0;

    wb1=waitbar(0,sprintf('EO tf continuation [%.2f %%] of %.0f/%.0f',0,id,length(TO_ref)));

    while ~iscomplete

        prob(L+it)=prob(id);

        if it==1        %-0NPCM--------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0
                
                prob(L+it).tf_ad=max(TO_ref(id).tf_ad,prob(id).tf_ad-(Dt_max/(2^(f-1)))*86400/TU);
                prob(L+it).tf=prob(L+it).tf_ad*TU+prob(L+it).t0;

                if rem(f,2)==1 && ~nang % attempt 0NPCM (alternates with act)

                    ll_g=prob(id).y0(8:14);        
                   
                else    % ACT

                    ll_g=ACT(prob(L+it));

                end

                if anynan(ll_g)
                    nang=1;
                    ex_flag=0;
                else
                    [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);       
                end

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
% 
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

                if rem(f,2)==1 && ~nang % attempt 1NPCM

                    ll_g=prob(L+it-1).y0(8:14)+(prob(L+it).tf_ad-prob(L+it-1).tf_ad)*(prob(L+it-1).y0(8:14)-prob(id).y0(8:14))/(prob(L+it-1).tf_ad-prob(id).tf_ad);            
                   
                else    % ACT

                    ll_g=ACT(prob(L+it));

                end

                if anynan(ll_g)
                    nang=1;
                    ex_flag=0;
                else
                    [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);       
                end                

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
                    [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(L+it)),ll_g,fsopt);       
                end                

                df=FO_ZFP(ll_FO,prob(L+it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                if prob(L+it).tf_ad-TO_ref(id).tf_ad<Dt_min*86400/TU && f>10
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

        nang=0;

        if prob(L+it).tf_ad==TO_ref(id).tf_ad

            iscomplete=1;
            prob(end).sts='TO';

        end

        if skip

            iscomplete=1;
            prob=prob(1:end-1);
            prob(end).sts='skp';

        end

        
        wb1=waitbar((prob(id).tf_ad-prob(L+it).tf_ad)/(prob(id).tf_ad-TO_ref(id).tf_ad),wb1,sprintf('EO tf continuation [%.2f %%] of %.0f/%.0f',(prob(id).tf_ad-prob(L+it).tf_ad)/(prob(id).tf_ad-TO_ref(id).tf_ad)*100,id,length(TO_ref)));

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

    fprintf('\n')

    close(wb1);

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

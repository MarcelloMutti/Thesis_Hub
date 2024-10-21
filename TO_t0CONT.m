function [prob] = TO_t0CONT(prob)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-8,'FunctionTolerance',1e-8,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    it=1;

    t_wo=prob(it).tw(1);
    t_wc=prob(it).tw(2);

    prob(it).t0=t_wo;

    Dt_max=5; % [days]
    Dt_iter=86400;
    N=25;

    while prob(it).t0<t_wc

        if it==1        %-1st-solution-through-T-continuation--------------
            
            cont='No';

            while strcmp(cont,'No')

                ex_flag=0;
                f=1;
                atmp=1; % number of rundown attempts
                tf_gv=2*pi*linspace(1,2,N); % rundown guesses

                wb1=waitbar(0,'Generating first TO solution');
    
                while ex_flag<=0
    
                    tf_g=tf_gv(f);
                    lltf_g=[ACT(prob(it)); tf_g];
    
                    [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,prob(it)),lltf_g,fsopt);

                    

                    wb1=waitbar(f/(atmp*N),wb1,'Generating first TO solution');

                    f=f+1;

                    if f>atmp*N
                        atmp=atmp+1;
                        tf_gv=2*pi*linspace(1,2,atmp*N);
                        f=1;
                    end
    
                end

                close(wb1);
                
                [~,~,prob(it)]=TO_ZFP(lltf_TO,prob(it));
                DispRes(prob(it),1);

                cont=questdlg('Accept initial solution?','Time cont','Yes','No','Exit','Exit');
                
                if strcmp(cont,'Exit')
                    close all
                    error('Continuation Aborted');
                elseif strcmp(cont,'No')
                    close all
                    clc
                end
            end

            close all
            clc

            wb2=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),'Initiating TO t0 continuation');

        elseif it==2    %-0NPCM--------------------------------------------
            
            tic

            ex_flag=0;
            f=1;

            while ex_flag<=0

                Dt_atmp=Dt_iter/(2^(f-1));
                prob(it).t0=prob(it-1).t0+Dt_atmp;

                if rem(f,2)==1 % attempt 0NPCM (alternates with act)

                    lltf_g=[prob(it-1).y0(8:14); prob(it-1).tf_ad];        
                   
                else    % ACT

                    lltf_g=[ACT(prob(it)); prob(it-1).tf_ad];

                end

                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,prob(it)),lltf_g,fsopt);

                

                df=TO_ZFP(lltf_TO,prob(it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end

                f=f+1;

            end


        else    %-1-3NPCM--------------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0

                Dt_atmp=Dt_iter/(2^(f-1));
                prob(it).t0=min(prob(it-1).t0+Dt_atmp,t_wc);

                if rem(f,4)==1 % attempt 1NPCM

                    lltf_g=[prob(it-1).y0(8:14); prob(it-1).tf_ad]+(prob(it).t0-prob(it-1).t0)*([prob(it-1).y0(8:14); prob(it-1).tf_ad]-[prob(it-2).y0(8:14); prob(it-2).tf_ad])/(prob(it-1).t0-prob(it-2).t0);            
                   
                elseif rem(f,4)==2 && it>=4 % attempt 2NPCM

                    yy=[prob(it-3:it-1).y0];
                    ll=yy(8:14,:);
                    lltf_g=makima([prob(it-3:it-1).t0],[ll; prob(it-3:it-1).tf_ad],prob(it).t0);

                elseif rem(f,4)==3 && it>=5 % attempt 3NPCM

                    yy=[prob(it-4:it-1).y0];
                    ll=yy(8:14,:);
                    lltf_g=makima([prob(it-4:it-1).t0],[ll; prob(it-4:it-1).tf_ad],prob(it).t0);

                else    % ACT

                    lltf_g=[ACT(prob(it)); prob(it-1).tf_ad];

                end

                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,prob(it)),lltf_g,fsopt);

                

                df=TO_ZFP(lltf_TO,prob(it));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end

                f=f+1;

            end

        end

        [~,~,prob(it)]=TO_ZFP(lltf_TO,prob(it));

        prob(it)=DispRes(prob(it),0);

        if prob(it).t0<t_wc
            prob(it+1)=prob(it);
            
            if it~=1
                Dt_iter=Dt_atmp;

                if (f-1)==1
                    Dt_iter=min(1.1*Dt_iter,Dt_max*86400);
                end
            end

            wb2=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),wb2,sprintf('TO t0 continuation [%.2f %%]',(prob(it).t0-t_wo)/(t_wc-t_wo)*100));

            it=it+1;
            
        end

        


        
    end

    fprintf('\n')

    close(wb2);

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
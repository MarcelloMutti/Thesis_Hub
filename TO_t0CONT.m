function [prob] = TO_t0CONT(prob)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'MaxIterations',2e2);

    it=1;

    t_wo=prob(it).tw(1);
    t_wc=prob(it).tw(2);

    prob(it).t0=t_wo;

    Dt=86400;
    N=25;

    while prob(it).t0<t_wc

        if it==1        %-1st-solution-through-T-continuation--------------
            
            cont='No';

            while strcmp(cont,'No')

                ex_flag=0;
                f=1;
                atmp=1; % number of rundown attempts
                tf_gv=linspace(2*pi,4*pi,N); % rundown guesses

                wb1=waitbar(0,'Generating first solution');
    
                while ex_flag<=0
    
                    tf_g=tf_gv(f);
                    lltf_g=[ACT(prob(it)); tf_g];
    
                    [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,prob(it)),lltf_g,fsopt);

                    wb1=waitbar(f/(atmp*N),wb1,'Generating first solution');

                    f=f+1;

                    if f>N
                        atmp=atmp+1;
                        tf_gv=linspace(2*pi,4*pi,atmp*N);
                        f=1;
                    end
    
                end

                close(wb1);
                
                [~,~,prob(it)]=TO_ZFP(lltf_TO,prob(it));
                DispRes(prob(it));

                cont=questdlg('Accept initial solution?','Time cont','Yes','No','Exit','Exit');
                
                if strcmp(cont,'Exit')
                    error('Continuation Aborted');
                end
            end

            close all
            clc

            wb2=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),'Initiating continuation');

        elseif it==2    %-0NPCM--------------------------------------------
            fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'MaxIterations',2e2);
            
            tic

            ex_flag=0;
            f=1;

            while ex_flag<=0

                Dt=Dt/(2^(f-1));
                prob(it).t0=prob(it-1).t0+Dt;

                if f==1 % attempt 0NPCM

                    lltf_g=[prob(it-1).y0(8:14); prob(it-1).tf_ad];        
                   
                else    % ACT

                    lltf_g=[ACT(prob(it)); prob(it-1).tf_ad];

                end

                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,prob(it)),lltf_g,fsopt);

                f=f+1;

            end


        else    %-1-3NPCM--------------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0

                Dt=Dt/(2^(f-1));
                prob(it).t0=min(prob(it-1).t0+Dt,t_wc);

                if f==1 % attempt 1NPCM

                    lltf_g=[prob(it-1).y0(8:14); prob(it-1).tf_ad]+(prob(it).t0-prob(it-1).t0)*([prob(it-1).y0(8:14); prob(it-1).tf_ad]-[prob(it-2).y0(8:14); prob(it-2).tf_ad])/(prob(it-1).t0-prob(it-2).t0);            
                   
                elseif f==2 && it>=4 % attempt 2NPCM

                    yy=[prob(it-3:it-1).y0];
                    ll=yy(8:14,:);
                    lltf_g=makima([prob(it-3:it-1).t0],[ll; prob(it-3:it-1).tf_ad],prob(it).t0);

                elseif f==3 && it>=5 % attempt 3NPCM

                    yy=[prob(it-4:it-1).y0];
                    ll=yy(8:14,:);
                    lltf_g=makima([prob(it-4:it-1).t0],[ll; prob(it-4:it-1).tf_ad],prob(it).t0);

                else    % ACT

                    lltf_g=[ACT(prob(it)); prob(it-1).tf_ad];

                end

                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,prob(it)),lltf_g,fsopt);

                f=f+1;

            end

%             fprintf('%.1f%%\n',(prob(it).t0-t_wo)/(t_wc-t_wo)*100)


        end

%         epsilon=0;

        [~,~,prob(it)]=TO_ZFP(lltf_TO,prob(it));
%         prob(it)=aux;

        prob(it)=DispRes(prob(it),0);

%         t0_v=[t0_v prob(it).t0];
%         lltf_M=[lltf_M lltf_TO];
%         lltf_alt=[lltf_alt [prob(it).y0(8:14); prob(it).tf_ad]];


        if prob(it).t0<t_wc
            prob(it+1)=prob(it);
            it=it+1;

            if f~=1 && it~=1    % check if correct
                Dt=min(1.1*Dt,5*86400);
            end
        end

        wb2=waitbar((prob(it).t0-t_wo)/(t_wc-t_wo),wb2,'TO continuation');


    end

    fprintf('\n')

    close(wb2);

    toc

end
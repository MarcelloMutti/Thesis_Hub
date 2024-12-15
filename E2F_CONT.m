function prob = E2F_CONT(prob,id,L)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    it=1;

    skip=0;
    nang=0;

    E_min=0;
    DE_max=0.05;
    DE=DE_max;

    DE_min=1e-6;
    E_skip=1e-5;
    % mp_tol=1e-8;

%     E=prob.epsilon;

    iscomplete=0;

    wb1=waitbar(0,sprintf('EO to FO continuation [%.2f %%] of %.0f/%.0f TO solution',0,id,L));

    while ~iscomplete

        prob(it+1)=prob(it);

        if it==1        %-0NPCM--------------------------------------------

            ex_flag=0;
            f=1;
            E=prob(it).epsilon;

            while ex_flag<=0

                prob(it+1).epsilon=max(E-DE/(2^(f-1)),E_min);

                if rem(f,2)==1 && ~nang % attempt 0NPCM (alternates with act)

                    ll_g=prob(it).y0(8:14);        
                   
                else    % ACT

                    ll_g=ACT(prob(it+1));

                end

                if anynan(ll_g)
                    nang=1;
                    ex_flag=0;
                else
                    [ll_FO,df,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it+1)),ll_g,fsopt);
                end                

                % df=FO_ZFP(ll_FO,prob(it+1));

                if (norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3) && prob(it+1).epsilon==0
                    ex_flag=0;
                end

                f=f+1;

            end


        else    %-1-3NPCM--------------------------------------------------

            ex_flag=0;
            f=1;
            E=prob(it).epsilon;

            while ex_flag<=0

                prob(it+1).epsilon=max(E-DE/(2^(f-1)),E_min);

                if rem(f,4+prob(it+1).isTO)==(1-(1-prob(it+1).isTO)) && ~nang % attempt 0NPCM

                    ll_g=prob(it).y0(8:14);

                elseif rem(f,4+prob(it+1).isTO)==(2-(1-prob(it+1).isTO)) && ~nang % attempt 1NPCM

                    ll_g=prob(it).y0(8:14)+(prob(it+1).epsilon-prob(it).epsilon)*(prob(it).y0(8:14)-prob(it-1).y0(8:14))/(prob(it).epsilon-prob(it-1).epsilon);            
                   
                elseif rem(f,4+prob(it+1).isTO)==(3-(1-prob(it+1).isTO)) && it>=3 && ~nang % attempt 2NPCM

                    yy=[prob(it-2:it).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(it-2:it).epsilon],ll,prob(it+1).epsilon);

                elseif rem(f,4+prob(it+1).isTO)==(4-(1-prob(it+1).isTO)) && it>=4 && ~nang % attempt 3NPCM

                    yy=[prob(it-3:it).y0];
                    ll=yy(8:14,:);
                    ll_g=makima([prob(it-3:it).epsilon],ll,prob(it+1).epsilon);

                else    % ACT

                    ll_g=ACT(prob(it+1));

                end

                if anynan(ll_g)
                    nang=1;
                    ex_flag=0;
                else
                    [ll_FO,df,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it+1)),ll_g,fsopt);
                end
                

                % df=FO_ZFP(ll_FO,prob(it+1));

                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3
                    ex_flag=0;
                end
                if (DE/(2^f)<DE_min && E<E_skip) || f>10
                    skip=1;
                    check=(ex_flag>0 && norm(df(1:3))*LU<10 && norm(df(4:6))*LU/TU<1e-3);
                    ex_flag=1;
                end

                f=f+1;

            end

        end


        [~,~,prob(it+1)]=FO_ZFP(ll_FO,prob(it+1));

        prob(it+1)=DispRes(prob(it+1),0);

        E=prob(it+1).epsilon;
        DE=max(min(E,DE_max),E_min);

        nang=0;

        % minimum step set to E_min, to check it ensures convergence

        if prob(it+1).epsilon==E_min || skip

            iscomplete=1;
            
        end

        if iscomplete==0
            wb1=waitbar((1-E)/(1-E_min),wb1,sprintf('EO to FO continuation [%.2f %%] of %.0f/%.0f TO solution',(1-E)/(1-E_min)*100,id,L));
        else
            wb1=waitbar(1,wb1,sprintf('EO to FO continuation [%.2f %%] of %.0f/%.0f TO solution',100,id,L));
        end

        it=it+1;

        


        
    end

    close(wb1);


    if skip==0
        prob=prob(end);
    else
        prob=prob(end-1+check);
    end


    
    

end

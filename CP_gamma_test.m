function [kk]=CP_gamma_test(prob,TO_ref,ep)

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9,'MaxIterations',2e2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    gL=zeros(size(TO_ref));
    kk=zeros(size(TO_ref));
    M=zeros(size(TO_ref));

    % S_Tol=1e-4;

    if ep==0
        wb1=waitbar(0,sprintf('FO(TO) testing of %.0f/%.0f',1,length(TO_ref)));
    else
        wb1=waitbar(0,sprintf('EO(TO) testing of %.0f/%.0f',1,length(TO_ref)));
    end

    for it=1:length(TO_ref)

        prob(it)=prob(end);

        %-analytical-rel-exploit---------------------------

        ex_flag=0;

        prob(it).tf_ad=TO_ref(it).tf_ad;
        prob(it).tf=TO_ref(it).tf;
        prob(it).t0=TO_ref(it).t0;
        prob(it).epsilon=ep;
        gL(it)=-1/max(TO_ref(it).S);

        k=1+ep;
        f=0;

        while ex_flag<=0

            ll_g=k*gL(it)*TO_ref(it).y0(8:14);

            [ll_FO,df,ex_flag]=fsolve(@(ll) FO_ZFP(ll,prob(it)),ll_g,fsopt);
            df_i=FO_ZFP(ll_g,prob(it));

            % df=FO_ZFP(ll_FO,prob(it));

            if k==ep+1
                if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3 || ex_flag<=0 || ~isequal(df,df_i)
                    k=ep+1.1;
                    ex_flag=0;
                else
                    ex_flag=1;
                end
            else
                if f==0
                    if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3 || ex_flag<=0
                        k=1.01*k;
                        ex_flag=0;
                    else
                        f=1;
                        k_old=k;
                        ll_old=ll_FO;
                        k=ep+1+(k-ep-1)*0.9;
                        ex_flag=0;
                    end
                else
                    if norm(df(1:3))*LU>10 || norm(df(4:6))*LU/TU>1e-3 || ex_flag<=0 || ~isequal(df,df_i) || k-(1+ep)<10*eps
                        k=k_old;
                        ll_FO=ll_old;
                        ex_flag=1;
                    else
                        k_old=k;
                        ll_old=ll_FO;
                        k=ep+1+(k-ep-1)*0.9;
                        ex_flag=0;
                    end
                end
            end

            % f=f+1;

        end
        
        if ep==0
            wb1=waitbar(it/length(TO_ref),wb1,sprintf('FO(TO) testing of %.0f/%.0f',min(it+1,length(TO_ref)),length(TO_ref)));
        else
            wb1=waitbar(it/length(TO_ref),wb1,sprintf('EO(TO) testing of %.0f/%.0f',min(it+1,length(TO_ref)),length(TO_ref)));
        end

        % g=mean(prob(it).y0(8:14)./TO_ref(it).y0(8:14));

        [~,~,prob(it)]=FO_ZFP(ll_FO,prob(it));
        prob(it)=DispRes(prob(it),0);

        kk(it)=k;
        M(it)=abs(ep+max(prob(it).S));
        
    end

    fprintf('\n')

    close(wb1);

    figure
    subplot(2,1,1)
    semilogy(et2MJD2000([TO_ref.t0]),kk,'-o')
    grid on
    grid minor
    axis tight
    ylabel('$$k=\frac{\gamma}{\gamma_L}\,[-]$$','interpreter','latex','fontsize',15)

    subplot(2,1,2)
    semilogy(et2MJD2000([TO_ref.t0]),M,'-o')
    grid on
    grid minor
    axis tight
    if ep==0
        ylabel('$$\|\max S_{FO}(t)-\varepsilon\|\,[-]$$','interpreter','latex','fontsize',15)
    else
        ylabel('$$\|\max S_{EO}(t)-\varepsilon\|\,[-]$$','interpreter','latex','fontsize',15)
    end
    xlabel('$$t_0\,[MJD2000]$$','interpreter','latex','fontsize',15)

end


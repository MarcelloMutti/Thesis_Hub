function lltf = TO_TCONT(t0,m0,targ)

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1    

    N=15;

    ToF_g0=100;
    M0=100;

    MaxIt=2e2;

%     fsopt=optimoptions('fsolve','Display','iter-detailed','FunctionTolerance',1e-12,'SpecifyObjectiveGradient',true);

    M=(M0-1)*(1-cos(linspace(pi/2,0,N)))+1;

%     sc_param_m=sc_param;

    lltf_M=zeros(8,N);

    wb=waitbar(0,'Intiating thrust continuation');

    for i=1:N

%         sc_param_m(1:2)=m(i)*sc_param(1:2);

        if i==1 % first solution attempt
            fsopt=optimoptions('fsolve','Display','iter-detailed','FunctionTolerance',1e-12,'SpecifyObjectiveGradient',true,'MaxIterations',MaxIt);


            ToF_g=ToF_g0;
            tf_g=ToF_g*86400/TU;

            ex_flag=0;
            while ex_flag<=0

%                 x0=SEL2_ND(t0);
%                 r0=norm(x0(1:3));
                
                Tc=M(i)*MARGO_param(1);

                l0_g=ACT(t0,[Tc; m0]);

                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,m0,targ,M(i)),[l0_g; tf_g],fsopt);
            end

            tf_g_old=tf_g;

        elseif i==2 % 0NPCM
            fsopt=optimoptions('fsolve','Display','iter-detailed','FunctionTolerance',1e-12,'SpecifyObjectiveGradient',true);

            ex_flag=0;
            f=0; % fail count
            tf_g=lltf_M(8,i-1);

            while ex_flag<=0

                if f==0 % attempt 0NPCM
                    l0_g=lltf_M(1:7,i-1);                
                else % ACT
                    l0_g=ACT(t0,sc_param_m);
                end

%                 fsopt=optimoptions('fsolve','Display','iter','FunctionTolerance',1e-12,'MaxFunctionEvaluations',(f+1)*maxFEval,'MaxIterations',(f+1)*maxIt);
                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,sc_param_m,targ),[l0_g; tf_g],fsopt);
                f=f+1;
            end

            tf_g_old=tf_g;

        elseif i==3 % 1NPCM
            ex_flag=0;
            f=0; % fail count
            tf_g=lltf_M(8,i-1)+(M(i)-M(i-1))*(lltf_M(8,i-1)-lltf_M(8,i-2))/(M(i-1)-M(i-2));

            while ex_flag<=0

                if f==0 % attempt linear
                    l0_g=lltf_M(1:7,i-1)+(M(i)-M(i-1))*(lltf_M(1:7,i-1)-lltf_M(1:7,i-2))/(M(i-1)-M(i-2));
                else % ACT
                    l0_g=ACT(t0,sc_param_m);
                end

%                 fsopt=optimoptions('fsolve','Display','iter','FunctionTolerance',1e-12,'MaxFunctionEvaluations',(f+1)*maxFEval,'MaxIterations',(f+1)*maxIt,'SpecifyObjectiveGradient',true);
                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,sc_param_m,targ),[l0_g; tf_g],fsopt);
                f=f+1;
            end

            tf_g_old=tf_g;

        else % power extrap on tf, up to 2NPCM on ll
            ex_flag=0;
            f=0; % fail count

%             tf_g=exp(polyval(polyfit(log(m(1:i-1)),log(lltf(8,1:i-1)),2),log(m(i))));
            tf_g=exp(polyval(polyfit(log(M(i-3:i-1)),log(lltf_M(8,i-3:i-1)),2),log(M(i))));

            tf_g=tf_g*lltf_M(8,i-1)/tf_g_old;

            while ex_flag<=0
                
                if f==0 % attempt linear
                    l0_g=lltf_M(1:7,i-1)+(M(i)-M(i-1))*(lltf_M(1:7,i-1)-lltf_M(1:7,i-2))/(M(i-1)-M(i-2));
                elseif f==1  % attempt makima 2nd order
                    l0_g=makima(M(i-3:i-1),lltf_M(1:7,i-3:i-1),M(i));
                else % ACT
                    l0_g=ACT(t0,sc_param_m);
                end

%                 fsopt=optimoptions('fsolve','Display','iter','FunctionTolerance',1e-12,'MaxFunctionEvaluations',(f+1)*maxFEval,'MaxIterations',(f+1)*maxIt,'SpecifyObjectiveGradient',true);
                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,sc_param_m,targ),[l0_g; tf_g],fsopt);
                f=f+1;
                
            end

            tf_g_old=tf_g;

        end


        lltf_M(:,i)=lltf_TO;


%         fprintf('%.1f%%\n',i/N*100)
        wb=waitbar(i/N,wb,'Thrust continuation');
    end
    fprintf('\n')
    close(wb);

    lltf=lltf_M(:,end);
end
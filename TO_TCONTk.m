function [m,lltf,dt] = TO_TCONTk(t0,sc_param,targ)

    N=25;
    ToF_g0=100;
    m0=100;

    maxFEval=1500;
    maxIt=500;

    fsopt=optimoptions('fsolve','Display','final','FunctionTolerance',1e-12,'MaxFunctionEvaluations',maxFEval,'MaxIterations',maxIt);
    m=(m0-1)*(1-cos(linspace(pi/2,0,N)))+1;
%     m=linspace(m0,1,N);
    sc_param_m=sc_param;

    lltf=zeros(8,N);
    dt=zeros(1,N);

    for i=1:N
        fprintf('%.0f/%.0f\n',i,N)
%         disp(i,'/',N)
        sc_param_m(1:2)=m(i)*sc_param(1:2);

        if i==1 % first solution attempt
            ToF_g=ToF_g0;
            tf_g=t0+ToF_g*86400;

            ex_flag=0;
            while ex_flag<=0
                l0_g=ACT(t0,sc_param_m);
                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,sc_param_m,targ),[l0_g; tf_g],fsopt);
            end
            dt(i)=abs(lltf_TO(8)-tf_g);
            tf_g_old=tf_g;

        elseif i==2 % 0NPCM
            ex_flag=0;
            f=0; % fail count
            tf_g=lltf(8,i-1);

            while ex_flag<=0

                if f==0 % attempt 0NPCM
                    l0_g=lltf(1:7,i-1);                
                else % ACT
                    l0_g=ACT(t0,sc_param_m);
                end

                fsopt=optimoptions('fsolve','Display','final','FunctionTolerance',1e-12,'MaxFunctionEvaluations',(f+1)*maxFEval,'MaxIterations',(f+1)*maxIt);
                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,sc_param_m,targ),[l0_g; tf_g],fsopt);
                f=f+1;
            end
            dt(i)=abs(lltf_TO(8)-tf_g);
            tf_g_old=tf_g;

        elseif i==3 % 1NPCM
            ex_flag=0;
            f=0; % fail count
            tf_g=lltf(8,i-1)+(m(i)-m(i-1))*(lltf(8,i-1)-lltf(8,i-2))/(m(i-1)-m(i-2));

            while ex_flag<=0

                if f==0 % attempt linear
                    l0_g=lltf(1:7,i-1)+(m(i)-m(i-1))*(lltf(1:7,i-1)-lltf(1:7,i-2))/(m(i-1)-m(i-2));
                else % ACT
                    l0_g=ACT(t0,sc_param_m);
                end

                fsopt=optimoptions('fsolve','Display','final','FunctionTolerance',1e-12,'MaxFunctionEvaluations',(f+1)*maxFEval,'MaxIterations',(f+1)*maxIt);
                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,sc_param_m,targ),[l0_g; tf_g],fsopt);
                f=f+1;
            end
            dt(i)=abs(lltf_TO(8)-tf_g);
            tf_g_old=tf_g;

        else % power extrap on tf, up to 2NPCM on ll
            ex_flag=0;
            f=0; % fail count

%             tf_g=exp(polyval(polyfit(log(m(1:i-1)),log(lltf(8,1:i-1)),2),log(m(i))));
            tf_g=exp(polyval(polyfit(log(m(i-3:i-1)),log(lltf(8,i-3:i-1)),1),log(m(i))));
            tf_g=tf_g*lltf(8,i-1)/tf_g_old;

            while ex_flag<=0

                
                if f==0 % attempt linear
                    l0_g=lltf(1:7,i-1)+(m(i)-m(i-1))*(lltf(1:7,i-1)-lltf(1:7,i-2))/(m(i-1)-m(i-2));
                elseif f==1  % attempt makima 2nd order
                    l0_g=makima(m(i-3:i-1),lltf(1:7,i-3:i-1),m(i));
                else % ACT
                    l0_g=ACT(t0,sc_param_m);
                end

                fsopt=optimoptions('fsolve','Display','final','FunctionTolerance',1e-12,'MaxFunctionEvaluations',(f+1)*maxFEval,'MaxIterations',(f+1)*maxIt);
                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,t0,sc_param_m,targ),[l0_g; tf_g],fsopt);
                f=f+1;
                
            end
            dt(i)=abs(lltf_TO(8)-tf_g);
            tf_g_old=tf_g;

        end


        lltf(:,i)=lltf_TO;

%         disp((lltf_TO(8)-t0)/86400)
% 
%         figure
%         plot(m(i),tf_g,'x')
%         hold on
%         plot(m(1:i),lltf(8,1:i),'-o')

    end
    dt=dt./86400;
end
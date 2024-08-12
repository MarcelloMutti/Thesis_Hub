function [t0_v,lltf_M,mp,Hf,mS] = TO_t0CONT(twin,sc_param,targ)

    fsopt=optimoptions('fsolve','Display','none','SpecifyObjectiveGradient',true);

    t_wo=twin(1);
    t_wc=twin(2);

    T0=t_wo;
    Dt=1*86400;

    it=1;

    t0_v=[];
    lltf_M=[];
    mp=[];
    Hf=[];
    mS=[];

    while T0<t_wc

        if it==1        %-1st-solution-through-T-continuation--------------

%             fprintf('Initiating first solution\n\n')

            lltf_TO=TO_TCONT(T0,sc_param,targ);

%             fprintf('Initiating continuation\n\n')
            wb=waitbar((T0-t_wo)/(t_wc-t_wo),'Initiating continuation');

        elseif it==2    %-0NPCM--------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0

                T0=t0_v(it-1)+Dt/f;

                if f==1 % attempt 0NPCM

                    lltf_g=lltf_M(:,it-1);            
                   
                else    % ACT

                    lltf_g=[ACT(T0,sc_param); lltf_M(8,it-1)];

                end

                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,T0,sc_param,targ),lltf_g,fsopt);

                f=f+1;

            end

%         elseif it==3    %-1NPCM--------------------------------------------
% 
%             ex_flag=0;
%             f=1;
% 
%             while ex_flag<=0
% 
%                 T0=t0_v(it-1)+Dt/f;
% 
%                 if f==1 % attempt 1NPCM
% 
%                     lltf_g=lltf_M(:,it-1)+(T0-t0_v(it-1))*(lltf_M(:,it-1)-lltf_M(:,it-2))/(t0_v(it-1)-t0_v(it-2));            
%                    
%                 else    % ACT
% 
%                     lltf_g=[ACT(T0,sc_param); lltf_M(8,it-1)];
% 
%                 end
% 
%                 [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,T0,sc_param,targ),lltf_g,fsopt);
% 
%                 f=f+1;
% 
%             end
% 
% 
%         elseif it==4    %-2NPCM--------------------------------------------
% 
%             ex_flag=0;
%             f=1;
% 
%             while ex_flag<=0
% 
%                 T0=t0_v(it-1)+Dt/f;
% 
%                 if f==1 % attempt 1NPCM
% 
%                     lltf_g=lltf_M(:,it-1)+(T0-t0_v(it-1))*(lltf_M(:,it-1)-lltf_M(:,it-2))/(t0_v(it-1)-t0_v(it-2));            
%                    
%                 elseif f==2 % attempt 2NPCM
% 
%                     lltf_g=makima(t0_v(it-1:it-3),lltf_M(:,it-1:it-3),T0);
% 
%                 else    % ACT
% 
%                     lltf_g=[ACT(T0,sc_param); lltf_M(8,it-1)];
% 
%                 end
% 
%                 [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,T0,sc_param,targ),lltf_g,fsopt);
% 
%                 f=f+1;
% 
%             end

        else    %-3NPCM----------------------------------------------------

            ex_flag=0;
            f=1;

            while ex_flag<=0

                T0=min(t0_v(it-1)+Dt/f,t_wc);

                if f==1 % attempt 1NPCM

                    lltf_g=lltf_M(:,it-1)+(T0-t0_v(it-1))*(lltf_M(:,it-1)-lltf_M(:,it-2))/(t0_v(it-1)-t0_v(it-2));            
                   
                elseif f==2 && it>=4 % attempt 2NPCM

                    lltf_g=makima(t0_v(it-3:it-1),lltf_M(:,it-3:it-1),T0);

                elseif f==3 && it>=5 % attempt 3NPCM

                    lltf_g=makima(t0_v(it-4:it-1),lltf_M(:,it-4:it-1),T0);

                else    % ACT

                    lltf_g=[ACT(T0,sc_param); lltf_M(8,it-1)];

                end

                [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,T0,sc_param,targ),lltf_g,fsopt);

                f=f+1;

            end

%             fprintf('%.1f%%\n',(T0-t_wo)/(t_wc-t_wo)*100)


        end

        epsilon=0;

        [~,yy,S,H]=DispRes(lltf_TO,T0,sc_param,targ,epsilon,0);

        t0_v=[t0_v T0];
        lltf_M=[lltf_M lltf_TO];

        mp=[mp (yy(1,7)-yy(end,7))*sc_param(3)];

        Hf=[Hf H(1)];

        mS=[mS max(S)];

        it=it+1;

        wb=waitbar((T0-t_wo)/(t_wc-t_wo),wb,'TO continuation');

    end

    fprintf('\n')

    close(wb);

end
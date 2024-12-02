function [tout,xout] = FO_ode87(prob,tspan,z0)

    
    % The coefficients of method
%     c=  [ 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1]';
%     
%     A = [ 1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
%               1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
%               1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
%               5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
%               3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0; 
%               29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
%               16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
%               39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
%               246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
%              -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
%               185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
%               403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0]';
%     
%      b8 = [ 14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731,   561292985/797845732,   -1041891430/1371343529,  760417239/1151165299, 118820643/751138087, -528747749/2220607170,  1/4]';
%     
%      b7 = [ 13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321, 656045339/265891186,   -3867574721/1518517206,   465885868/322736535,  53011238/667516719,                  2/45,    0]';
    

    fsopt=optimoptions("fsolve",'SpecifyObjectiveGradient',true,'FunctionTolerance',eps,'OptimalityTolerance',eps,'Display','none');
    fzopt=optimset('Display','notify','FunValCheck','on');
    
    pow = 1/8; % power for step control
    
    % Check inputs
    if nargin < 3
        error('Not enough input arguments')
    end
    
    % Maximal step size
    hmax = (tspan(2) - tspan(1))/100;
    
    % Minimal step size
    hmin = 16*eps;
    
    % initial step size
    h = hmax;
    % if h>0.1
    %   h=0.1;
    % elseif h>hmax 
    %   h = hmax;
    % end
    
    %  A relative error tolerance that applies to all components of the solution vector. 
    Tol=1e-12;
    
    t0 = tspan(1);
    tf = tspan(2);
    
    % constant for step rejection
    reject = 0;
    
    t = t0;
    z = z0(:);          % start point
    % f = z*zeros(1,13);  % array f for RHS calculation
    tout = t;
    xout = z.';
    
    % initial sw eval-----------------------
    r=norm(z(1:3));
    [~,~,Sp]=MARGO_param(r);
    if Sp<prob.Plim(2)
        Ptype='med';
    else
        Ptype='max';
    end

    ep=prob.epsilon;

    Se=SwFun(t,z,prob.isFO);
    if ep~=0
        if Se<-ep
            utype='on';
        elseif Se>ep
            utype='off';
        else
            utype='med';
        end
    else
        if Se<0
            utype='on';
        else
            utype='off';
        end
    end
    
    Ptype_old=Ptype;
    utype_old=utype;
    
    % The main loop
    
    while (t < tf) && (h >= hmin)
    
        if (t + h) > tf 
            h = tf - t; 
        end
    
    
%         % Compute the RHS for step of method
%         f(:,1) = TO_2BP_SEP(t,z,Ptype_old);
%         for j = 1: 12
%             f(:,j+1) = TO_2BP_SEP(t+c(j)*h, z+h*f*A(:,j),Ptype_old);
%         end
%     
%         % Two solution 
%         x8=z+h*f*b8;
%         x7=z+h*f*b7;

        [x7,x8]=step(t,z,h,Ptype_old,utype_old,ep);

        % Truncation error 
        err = norm(x7-x8);
      
        % Estimate the error and the acceptable error
        step_err = norm(err,'inf');
        tau = Tol*max(norm(z,'inf'),1.0);
    
        % Update the solution only if the error is acceptable
        if step_err <= tau

            % step sw eval @ t+h ------------------------------------------
            r=norm(x8(1:3));
            [~,~,Sp]=MARGO_param(r);
            if Sp<prob.Plim(2)
                Ptype='med';
            else
                Ptype='max';
            end

            r_old=norm(z(1:3));
            [~,~,Sp_old]=MARGO_param(r_old);

            % c=Tc(2);
            % cp=Tcp(2);
            % dtSe=c/x8(7)*dot(x8(8:10),x8(11:13))/norm(x8(11:13))-norm(x8(11:13))/x8(7)*cp/r*dot(x8(1:3),x8(4:6));

            Se_old=SwFun(t,z,prob.isFO);
            Se=SwFun(t+h,x8,prob.isFO);
            if ep~=0
                if Se<-ep
                    utype='on';
                elseif Se>ep
                    utype='off';
                else
                    utype='med';
                end
            else
                if Se<0
                    utype='on';
                else
                    utype='off';
                end
            end

            if (~strcmp(Ptype,Ptype_old) && ~strcmp(utype,utype_old)) ...   % Double switch
                || (ep~=0 && strcmp(utype,'on') && strcmp(utype_old,'off')) ... % off->on in EO
                || (ep~=0 && strcmp(utype,'off') && strcmp(utype_old,'on')) ... % on->off in EO
                || ((~strcmp(Ptype,Ptype_old) || ~strcmp(utype,utype_old)) && (isapprox(abs(Se_old),ep,'tight') || isapprox(Sp_old,prob.Plim(2),'tight')))  % force 1step after switching
                % || isequal(abs(Se),ep) ...                            % Sw=+-ep
                % || (isapprox(abs(dtSe),0,'tight') && isapprox(abs(Se),ep,'tight')) ... % dubious utility
                

                h=h*0.75;

                crossing=1;

            elseif xor(~strcmp(Ptype,Ptype_old),~strcmp(utype,utype_old))

                % bisection attempt (fzero)
                ex_flag=0;

                while ex_flag<=0

                    tc=[];

                    if ~strcmp(Ptype,Ptype_old) && power_switching(t,t,z,Ptype_old,utype_old,prob)*power_switching(t+h,t,z,Ptype_old,utype_old,prob)<0

                        [tc,cr_v,ex_flag]=fzero(@(T) power_switching(T,t,z,Ptype_old,utype_old,prob),[t, t+h],fzopt);

                    elseif ~strcmp(utype,utype_old) && throttle_switching(t,t,z,Ptype_old,utype_old,utype,prob)*throttle_switching(t+h,t,z,Ptype_old,utype_old,utype,prob)<0

                        [tc,cr_v,ex_flag]=fzero(@(T) throttle_switching(T,t,z,Ptype_old,utype_old,utype,prob),[t, t+h],fzopt);

                    end

                    if isempty(tc) || (tc<t || tc>t+h)
                        
                        % newton attempt (fsolve)
                        ex_flag=0;
                        tg=t+h/2;
        
                        while ex_flag<=0
        
                            if ~strcmp(Ptype,Ptype_old)
        
                                [tc,cr_v,ex_flag]=fsolve(@(T) power_switching(T,t,z,Ptype_old,utype_old,prob),tg,fsopt);
        
                            elseif ~strcmp(utype,utype_old)
        
                                [tc,cr_v,ex_flag]=fsolve(@(T) throttle_switching(T,t,z,Ptype_old,utype_old,utype,prob),tg,fsopt);
        
                            end
        
                            if tc<t || tc>t+h

                                % error('fsolve fail')
                                ex_flag=0;
                                tg=unifrnd(t,t+h);
        
                            end
        
                        end

                    end

                end

                %----begin switch block

                [~,zc]=step(t,z,tc-t,Ptype_old,utype_old,ep);           % step to crossing

                % perform tentative step from zc tc

                [~,z_pred]=step(tc,zc,t+h-tc,Ptype_old,utype_old,ep);   % predictive step

                [~,~,Sp_p]=MARGO_param(norm(z_pred(1:3)));
                if Sp_p<prob.Plim(2)
                    Ptype_p='med';
                else
                    Ptype_p='max';
                end

                Se_p=SwFun(t+h,z_pred,prob.isFO);
                if ep~=0
                    if Se_p<-ep
                        utype_p='on';
                    elseif Se_p>ep
                        utype_p='off';
                    else
                        utype_p='med';
                    end
                else
                    if Se_p<0
                        utype_p='on';
                    else
                        utype_p='off';
                    end
                end

                if (~strcmp(utype,utype_old) && strcmp(utype_p,utype_old)) || (~strcmp(Ptype,Ptype_old) && strcmp(Ptype_p,Ptype_old))
                    % tangency without switch

                    tout=[tout; tc; t+h];
                    xout=[xout; zc.'; z_pred.'];

                    t=t+h;
                    z=z_pred;

                else
    
                    rc=zc(1:3);
                    vc=zc(4:6);
                    mc=zc(7);
                    llrc=zc(8:10);
                    llvc=zc(11:13);
        
                    dzc_m=FO_2BP_SEP(tc,zc,Ptype_old,utype_old,ep);
                    dzc_p=FO_2BP_SEP(tc,zc,Ptype,utype,ep);
    
                    if ~strcmp(Ptype,Ptype_old)
        
                        Psi=eye(14)+(dzc_p(1:14)-dzc_m(1:14))*[rc.'/dot(rc,vc),zeros(1,11)];
    
                    elseif (~strcmp(utype,utype_old) && ep==0)
    
                        [Tc,Tcp]=MARGO_param(norm(rc));
                        c=Tc(2);    % ex vel
                        cp=Tcp(2);  % dc/dr
    
                        DySe=[-norm(llvc)/mc*cp*rc.'/norm(rc), zeros(1,3), c*norm(llvc)/mc^2, zeros(1,3), -c/mc*llvc.'/norm(llvc), -1];
                        DtSe=c/mc*dot(llrc,llvc)/norm(llvc)-norm(llvc)*cp/mc*dot(rc,vc)/norm(rc);
    
                        Psi=eye(14)+(dzc_p(1:14)-dzc_m(1:14))*DySe/DtSe;
    
                    else
    
                        Psi=eye(14);
    
                    end
        
                    h=tc-t;
    
                    % if abs(tc-t)<Tol
                    %     tang=1;
                    % else
                    %     tang=0;
                    % end
    
                    t=tc;
                    z=zc;
                    zp=z;
        
                    Phi_m=reshape(z(15:210),[14,14]);
                    Phi_p=Psi*Phi_m;
                    zp(15:210)=reshape(Phi_p,[14*14,1]);
        
    
                    % if tang
                    % 
                    %     [~,x8]=step(tout(end),xout(end-1,:).',h,Ptype,utype,ep);
                    % 
                    %     tout = [tout; t; t; tout(end)+h];
                    %     xout = [xout; z.'; zp.'; x8.'];
                    % 
                    %     t=tout(end);
                    %     z=x8;
                    % 
                    % else
    
                        tout = [tout; t; t];
                        xout = [xout; z.'; zp.'];
        
                        z=zp;
    
                    % end

                    Ptype_old=Ptype;
                    utype_old=utype;

                end
    
                crossing=1;

                %----end switch block
    
            else    % no switching detected
        
                t = t + h;
                z = x8; 
                tout = [tout; t];
                xout = [xout; z.'];
        
                reject = 0;

                crossing = 0;

            end
    
        else
    
            crossing=0;
            reject = 1;

        end
    
        
    
        % Step control
        if step_err == 0.0
            step_err = eps*10.0;
        end

        if crossing==0 && ~isnan(step_err)
            h = min(hmax, 0.9*h*(tau/step_err)^pow);
        elseif crossing==0 && isnan(step_err)
            h = h*0.75;
        end

        if (abs(h) <= eps) 
            if reject == 0
                disp('Warning!!! ode87. Step is very small!!!');
                h = eps * 100;
                return
            else
                disp('Error in ode87. Step is too small.');
                return 
            end
        end
    
    end
    
    if (t < tf)
        disp('Error in ODE87...')
        t
    end

end

function [x7,x8]=step(t,z,h,Ptype,utype,ep)

    C=  [ 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1]';
    
    A = [ 1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
              1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
              1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
              5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
              3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0; 
              29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
              16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
              39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
              246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
             -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
              185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
              403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0]';
    
     b8 = [ 14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731,   561292985/797845732,   -1041891430/1371343529,  760417239/1151165299, 118820643/751138087, -528747749/2220607170,  1/4]';
    
     b7 = [ 13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321, 656045339/265891186,   -3867574721/1518517206,   465885868/322736535,  53011238/667516719,                  2/45,    0]';


    f = z*zeros(1,13);
    % Compute the RHS for step of method
    f(:,1) = FO_2BP_SEP(t,z,Ptype,utype,ep);
    for j = 1:12
        f(:,j+1) = FO_2BP_SEP(t+C(j)*h, z+h*f*A(:,j),Ptype,utype,ep);
    end

    % Two solution 
    x8=z+h*f*b8;
    x7=z+h*f*b7;

end


function [f,df]=power_switching(tc,tk,zk,Ptype,utype,prob)

    h=tc-tk;
    
    [~,zc]=step(tk,zk,h,Ptype,utype,prob.epsilon);

    rc=zc(1:3);
    vc=zc(4:6);

    [~,~,Sp,dPdr]=MARGO_param(norm(rc));

    df=dPdr*dot(rc,vc)/norm(rc);

    f=Sp-prob.Plim(2);

end

function [f,df]=throttle_switching(tc,tk,zk,Ptype,utype,utype_new,prob)

    h=tc-tk;
    
    [~,zc]=step(tk,zk,h,Ptype,utype,prob.epsilon);

    rc=zc(1:3);
    vc=zc(4:6);
    mc=zc(7);
    llrc=zc(8:10);
    llvc=zc(11:13);

    [Tc,Tcp]=MARGO_param(norm(rc));
    c=Tc(2);    % ex vel
    cp=Tcp(2);  % dc/dr

    Se=SwFun(tc,zc,prob.isFO);

    ep=prob.epsilon;

    if ep==0    % FO case

        f=Se;

    else

        if (strcmp(utype,'on') && strcmp(utype_new,'med')) ...       % on->med in EO
                || (strcmp(utype,'med') && strcmp(utype_new,'on'))   % med->on in EO

            f=Se+ep;

        else    % off->med or med->off in EO

            f=Se-ep;

        end

    end

    df=c/mc*dot(llrc,llvc)/norm(llvc)-norm(llvc)*cp/mc*dot(rc,vc)/norm(rc);

end




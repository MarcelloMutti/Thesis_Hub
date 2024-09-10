function prob = struct_assembly(targ,twin,m0,Plim,epsilon)

    % targ    [str]
    % m0      [1x1, kg]
    % Plim    [1x2, W, (Pmin Pmax)]
    % twin    [1x2, s, (t_wo t_wc)]
    % epsilon [1x1]
    
    %-problem-setup--------------------------------------------------------
    prob.targ=targ;         % [str]
    prob.tw=twin;           % [1x2, s]
    prob.m0=m0;             % [1x1, kg]
    prob.Plim=Plim;         % [1x2, W]
    prob.epsilon=epsilon;   % [1x1, -]
    
    %-added-by-TO_t0CONT---------------------------------------------------
    prob.t0=[];             % [1x1, s]
    
    %-added-by-TO_ZFP------------------------------------------------------
    prob.gamma=[];          % [8x1, -]
    prob.jac=[];            % [8x8, -]
    prob.tf=[];             % [1x1, s]
    prob.tf_ad=[];          % [1x1, -]
    prob.y0=[];             % [14x1, -]

    %-added-by-DispRes-----------------------------------------------------
    prob.zz=[];             % [Nx210, -]
    prob.tt=[];             % [Nx1, d]
    prob.tt_ad=[];          % [Nx1, -]
    prob.mf=[];             % [1x1, kg]
    prob.mp=[];             % [1x1, kg]
    prob.S=[];              % [Nx1, -]
    prob.H=[];              % [Nx1, -]
    
end
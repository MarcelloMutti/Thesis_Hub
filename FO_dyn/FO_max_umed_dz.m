function [FF_umed,A_umed] = FO_max_umed_dz(in1,in2,in3,g0,epsilon)
%FO_max_umed_dz
%    [FF_umed,A_umed] = FO_max_umed_dz(IN1,IN2,IN3,G0,EPSILON)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    02-Dec-2024 17:48:31

MP1_1 = in2(1);
MP1_2 = in2(4);
MP1_3 = in2(7);
MP1_4 = in2(10);
MP1_5 = in2(13);
MP2_1 = in2(2);
MP2_2 = in2(5);
MP2_3 = in2(8);
MP2_4 = in2(11);
MP2_5 = in2(14);
U1 = in3(:,1);
U2 = in3(:,2);
U3 = in3(:,3);
U4 = in3(:,4);
y1 = in1(1,:);
y2 = in1(2,:);
y3 = in1(3,:);
y4 = in1(4,:);
y5 = in1(5,:);
y6 = in1(6,:);
y7 = in1(7,:);
y8 = in1(8,:);
y9 = in1(9,:);
y10 = in1(10,:);
y11 = in1(11,:);
y12 = in1(12,:);
y13 = in1(13,:);
y14 = in1(14,:);
t2 = abs(y1);
t3 = abs(y2);
t4 = abs(y3);
t5 = abs(y11);
t6 = abs(y12);
t7 = abs(y13);
t8 = sign(y1);
t9 = sign(y2);
t10 = sign(y3);
t11 = sign(y11);
t12 = sign(y12);
t13 = sign(y13);
t14 = U2.^2;
t15 = U2.^3;
t22 = y1.*y11.*3.0;
t23 = y2.*y12.*3.0;
t24 = y3.*y13.*3.0;
t25 = 1.0./U1;
t27 = 1.0./U3;
t28 = 1.0./epsilon;
t29 = 1.0./g0;
t30 = 1.0./y7;
t33 = MP1_2.*1.2e+2;
t34 = MP2_2.*1.2e+2;
t35 = MP1_3.*1.44e+4;
t36 = MP2_3.*1.44e+4;
t37 = MP1_4.*1.728e+6;
t38 = MP2_4.*1.728e+6;
t39 = MP1_5.*2.0736e+8;
t40 = MP2_5.*2.0736e+8;
t16 = t2.^2;
t17 = t3.^2;
t18 = t4.^2;
t19 = t5.^2;
t20 = t6.^2;
t21 = t7.^2;
t26 = t25.^2;
t31 = t30.^2;
t32 = t30.^3;
t43 = t22+t23+t24;
t60 = MP1_1+t33+t35+t37+t39;
t61 = MP2_1+t34+t36+t38+t40;
t41 = t16+t17+t18;
t42 = t19+t20+t21;
t62 = 1.0./t61;
t44 = 1.0./t42;
t45 = sqrt(t42);
t46 = 1.0./t41.^(3.0./2.0);
t47 = 1.0./t41.^(5.0./2.0);
t48 = 1.0./t41.^(7.0./2.0);
t49 = 1.0./t45;
t51 = -t46;
t52 = t47.*y1.*y2.*3.0;
t53 = t47.*y1.*y3.*3.0;
t54 = t47.*y2.*y3.*3.0;
t58 = t43.*t47;
t63 = U2.*g0.*t25.*t30.*t45.*t61;
t65 = (U4.*t14.*t25.*t27.*t28.*t31.*t45.*t60)./2.0;
t50 = t49.^3;
t55 = -t52;
t56 = -t53;
t57 = -t54;
t59 = -t58;
t64 = epsilon+t63+y14-1.0;
FF_umed = [y4;y5;y6;t51.*y1-(U4.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*t64.*y11)./2.0;t51.*y2-(U4.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*t64.*y12)./2.0;t51.*y3-(U4.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*t64.*y13)./2.0;U2.*U4.*t27.*t28.*t29.*t60.*t62.*t64.*(-1.0./2.0);t46.*y11+t59.*y1;t46.*y12+t59.*y2;t46.*y13+t59.*y3;-y8;-y9;-y10;U4.*t14.*t25.*t27.*t28.*t31.*t45.*t60.*t64.*(-1.0./2.0)];
if nargout > 1
    t66 = (U4.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*t64)./2.0;
    t67 = -t66;
    mt1 = [0.0,0.0,0.0,t51+t2.*t8.*t47.*y1.*3.0,t2.*t8.*t47.*y2.*3.0,t2.*t8.*t47.*y3.*3.0,0.0,t59-t47.*y1.*y11.*3.0-t2.*t8.*t47.*y11.*3.0+t2.*t8.*t43.*t48.*y1.*5.0,t47.*y2.*y11.*-3.0-t2.*t8.*t47.*y12.*3.0+t2.*t8.*t43.*t48.*y2.*5.0,t47.*y3.*y11.*-3.0-t2.*t8.*t47.*y13.*3.0+t2.*t8.*t43.*t48.*y3.*5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3.*t9.*t47.*y1.*3.0,t51+t3.*t9.*t47.*y2.*3.0,t3.*t9.*t47.*y3.*3.0,0.0,t47.*y1.*y12.*-3.0-t3.*t9.*t47.*y11.*3.0+t3.*t9.*t43.*t48.*y1.*5.0,t59-t47.*y2.*y12.*3.0-t3.*t9.*t47.*y12.*3.0+t3.*t9.*t43.*t48.*y2.*5.0,t47.*y3.*y12.*-3.0-t3.*t9.*t47.*y13.*3.0+t3.*t9.*t43.*t48.*y3.*5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t10.*t47.*y1.*3.0,t4.*t10.*t47.*y2.*3.0];
    mt2 = [t51+t4.*t10.*t47.*y3.*3.0,0.0,t47.*y1.*y13.*-3.0-t4.*t10.*t47.*y11.*3.0+t4.*t10.*t43.*t48.*y1.*5.0,t47.*y2.*y13.*-3.0-t4.*t10.*t47.*y12.*3.0+t4.*t10.*t43.*t48.*y2.*5.0,t59-t47.*y3.*y13.*3.0-t4.*t10.*t47.*y13.*3.0+t4.*t10.*t43.*t48.*y3.*5.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,(U4.*g0.*t15.*t26.*t27.*t28.*t32.*t60.*t61.*y11)./2.0+(U4.*t14.*t25.*t27.*t28.*t31.*t49.*t60.*t64.*y11)./2.0,(U4.*g0.*t15.*t26.*t27.*t28.*t32.*t60.*t61.*y12)./2.0+(U4.*t14.*t25.*t27.*t28.*t31.*t49.*t60.*t64.*y12)./2.0];
    mt3 = [(U4.*g0.*t15.*t26.*t27.*t28.*t32.*t60.*t61.*y13)./2.0+(U4.*t14.*t25.*t27.*t28.*t31.*t49.*t60.*t64.*y13)./2.0,t65,0.0,0.0,0.0,0.0,0.0,0.0,U4.*t14.*t25.*t27.*t28.*t32.*t45.*t60.*t64+(U4.*g0.*t15.*t26.*t27.*t28.*t31.^2.*t42.*t60.*t61)./2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,t67+(U4.*t5.*t11.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y11)./2.0-(U4.*g0.*t5.*t11.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y11)./2.0,(U4.*t5.*t11.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y12)./2.0-(U4.*g0.*t5.*t11.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y12)./2.0];
    mt4 = [(U4.*t5.*t11.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y13)./2.0-(U4.*g0.*t5.*t11.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y13)./2.0,U4.*t5.*t11.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*(-1.0./2.0),t46-t47.*y1.^2.*3.0,t55,t56,0.0,0.0,0.0,U4.*g0.*t5.*t11.*t15.*t26.*t27.*t28.*t32.*t60.*t61.*(-1.0./2.0)-(U4.*t5.*t11.*t14.*t25.*t27.*t28.*t31.*t49.*t60.*t64)./2.0,0.0,0.0,0.0,(U4.*t6.*t12.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y11)./2.0-(U4.*g0.*t6.*t12.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y11)./2.0,t67+(U4.*t6.*t12.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y12)./2.0-(U4.*g0.*t6.*t12.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y12)./2.0];
    mt5 = [(U4.*t6.*t12.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y13)./2.0-(U4.*g0.*t6.*t12.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y13)./2.0,U4.*t6.*t12.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*(-1.0./2.0),t55,t46-t47.*y2.^2.*3.0,t57,0.0,0.0,0.0,U4.*g0.*t6.*t12.*t15.*t26.*t27.*t28.*t32.*t60.*t61.*(-1.0./2.0)-(U4.*t6.*t12.*t14.*t25.*t27.*t28.*t31.*t49.*t60.*t64)./2.0,0.0,0.0,0.0,(U4.*t7.*t13.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y11)./2.0-(U4.*g0.*t7.*t13.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y11)./2.0,(U4.*t7.*t13.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y12)./2.0-(U4.*g0.*t7.*t13.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y12)./2.0];
    mt6 = [t67+(U4.*t7.*t13.*t14.*t25.*t27.*t28.*t30.*t50.*t60.*t64.*y13)./2.0-(U4.*g0.*t7.*t13.*t15.*t26.*t27.*t28.*t31.*t44.*t60.*t61.*y13)./2.0,U4.*t7.*t13.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*(-1.0./2.0),t56,t57,t46-t47.*y3.^2.*3.0,0.0,0.0,0.0,U4.*g0.*t7.*t13.*t15.*t26.*t27.*t28.*t32.*t60.*t61.*(-1.0./2.0)-(U4.*t7.*t13.*t14.*t25.*t27.*t28.*t31.*t49.*t60.*t64)./2.0,0.0,0.0,0.0,U4.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*y11.*(-1.0./2.0),U4.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*y12.*(-1.0./2.0),U4.*t14.*t25.*t27.*t28.*t30.*t49.*t60.*y13.*(-1.0./2.0),U2.*U4.*t27.*t28.*t29.*t60.*t62.*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,-t65];
    A_umed = reshape([mt1,mt2,mt3,mt4,mt5,mt6],14,14);
end
end

function [FF_med_uon,A_med_uon] = FO_med_uon_dz(in1,in2,in3,g0,epsilon)
%FO_med_uon_dz
%    [FF_med_uon,A_med_uon] = FO_med_uon_dz(IN1,IN2,IN3,G0,EPSILON)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    14-Nov-2024 00:25:05

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
MP3_1 = in2(3);
MP3_2 = in2(6);
MP3_3 = in2(9);
MP3_4 = in2(12);
MP3_5 = in2(15);
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
t21 = y1.*y11.*3.0;
t22 = y2.*y12.*3.0;
t23 = y3.*y13.*3.0;
t24 = 1.0./U1;
t25 = 1.0./U2;
t26 = 1.0./U3;
t27 = 1.0./g0;
t28 = y14-1.0;
t29 = 1.0./y7;
t15 = t2.^2;
t16 = t3.^2;
t17 = t4.^2;
t18 = t5.^2;
t19 = t6.^2;
t20 = t7.^2;
t30 = t29.^2;
t31 = MP3_3.*t2.*t8.*2.0;
t32 = MP3_3.*t3.*t9.*2.0;
t33 = MP3_3.*t4.*t10.*2.0;
t34 = MP3_4.*t2.*t8.*6.0;
t35 = MP3_4.*t3.*t9.*6.0;
t36 = MP3_4.*t4.*t10.*6.0;
t39 = t21+t22+t23;
t37 = t15+t16+t17;
t38 = t18+t19+t20;
t40 = t37.^2;
t41 = MP3_3.*t37;
t42 = MP3_4.*t37.*3.0;
t44 = sqrt(t37);
t46 = sqrt(t38);
t57 = MP3_5.*t2.*t8.*t37.*4.0;
t58 = MP3_5.*t3.*t9.*t37.*4.0;
t59 = MP3_5.*t4.*t10.*t37.*4.0;
t43 = MP3_5.*t40;
t45 = t44.^3;
t47 = 1.0./t44;
t51 = 1.0./t46;
t53 = MP3_2.*t44;
t55 = MP3_3.*t44.*2.0;
t70 = MP3_4.*t2.*t8.*t44.*3.0;
t71 = MP3_4.*t3.*t9.*t44.*3.0;
t72 = MP3_4.*t4.*t10.*t44.*3.0;
t76 = MP3_5.*t2.*t8.*t44.*1.2e+1;
t77 = MP3_5.*t3.*t9.*t44.*1.2e+1;
t78 = MP3_5.*t4.*t10.*t44.*1.2e+1;
t48 = 1.0./t45;
t49 = t47.^5;
t50 = t47.^7;
t52 = t51.^3;
t54 = MP3_4.*t45;
t56 = MP3_5.*t45.*4.0;
t67 = MP3_2.*t2.*t8.*t47;
t68 = MP3_2.*t3.*t9.*t47;
t69 = MP3_2.*t4.*t10.*t47;
t73 = t31.*t47;
t74 = t32.*t47;
t75 = t33.*t47;
t60 = -t48;
t61 = t49.*y1.*y2.*3.0;
t62 = t49.*y1.*y3.*3.0;
t63 = t49.*y2.*y3.*3.0;
t79 = t39.*t49;
t81 = t34+t73+t76;
t82 = t35+t74+t77;
t83 = t36+t75+t78;
t84 = MP3_2+t42+t55+t56;
t85 = MP3_1+t41+t43+t53+t54;
t93 = t31+t57+t67+t70;
t94 = t32+t58+t68+t71;
t95 = t33+t59+t69+t72;
t64 = -t61;
t65 = -t62;
t66 = -t63;
t80 = -t79;
t86 = t85.^2;
t87 = t85.^3;
t89 = MP1_2.*t85;
t90 = MP2_2.*t85;
t91 = MP1_3.*t85.*2.0;
t92 = MP2_3.*t85.*2.0;
t106 = MP1_2.*t93;
t107 = MP1_2.*t94;
t108 = MP2_2.*t93;
t109 = MP1_2.*t95;
t110 = MP2_2.*t94;
t111 = MP2_2.*t95;
t112 = MP1_3.*t93.*2.0;
t113 = MP1_3.*t94.*2.0;
t114 = MP2_3.*t93.*2.0;
t115 = MP1_3.*t95.*2.0;
t116 = MP2_3.*t94.*2.0;
t117 = MP2_3.*t95.*2.0;
t119 = MP1_4.*t85.*t93.*6.0;
t122 = MP1_4.*t85.*t94.*6.0;
t123 = MP2_4.*t85.*t93.*6.0;
t126 = MP1_4.*t85.*t95.*6.0;
t127 = MP2_4.*t85.*t94.*6.0;
t129 = MP2_4.*t85.*t95.*6.0;
t88 = t86.^2;
t96 = MP1_3.*t86;
t97 = MP1_4.*t87;
t99 = MP2_3.*t86;
t100 = MP2_4.*t87;
t102 = MP1_4.*t86.*3.0;
t103 = MP1_5.*t87.*4.0;
t104 = MP2_4.*t86.*3.0;
t105 = MP2_5.*t87.*4.0;
t118 = t91.*t93;
t120 = t91.*t94;
t121 = t92.*t93;
t124 = t91.*t95;
t125 = t92.*t94;
t128 = t92.*t95;
t142 = MP1_5.*t86.*t93.*1.2e+1;
t143 = MP1_5.*t86.*t94.*1.2e+1;
t144 = MP2_5.*t86.*t93.*1.2e+1;
t145 = MP1_5.*t86.*t95.*1.2e+1;
t146 = MP2_5.*t86.*t94.*1.2e+1;
t147 = MP2_5.*t86.*t95.*1.2e+1;
t98 = MP1_5.*t88;
t101 = MP2_5.*t88;
t130 = t93.*t102;
t131 = t93.*t103;
t132 = t94.*t102;
t133 = t93.*t104;
t134 = t94.*t103;
t135 = t93.*t105;
t136 = t95.*t102;
t137 = t94.*t104;
t138 = t95.*t103;
t139 = t94.*t105;
t140 = t95.*t104;
t141 = t95.*t105;
t148 = MP1_2+t91+t102+t103;
t149 = MP2_2+t92+t104+t105;
t161 = t112+t119+t142;
t162 = t113+t122+t143;
t163 = t114+t123+t144;
t164 = t115+t126+t145;
t165 = t116+t127+t146;
t166 = t117+t129+t147;
t150 = MP1_1+t89+t96+t97+t98;
t151 = MP2_1+t90+t99+t100+t101;
t154 = U4.*t14.*t24.*t26.*t47.*t84.*t148;
t167 = t106+t118+t130+t131;
t168 = t107+t120+t132+t134;
t169 = t108+t121+t133+t135;
t170 = t109+t124+t136+t138;
t171 = t110+t125+t137+t139;
t172 = t111+t128+t140+t141;
t152 = 1.0./t151;
t155 = t154.*y1;
t156 = t154.*y2;
t157 = t154.*y3;
t158 = U4.*t14.*t24.*t26.*t29.*t51.*t150;
t160 = t29.*t46.*t154;
t153 = t152.^2;
t159 = -t158;
t173 = U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t152;
t174 = t173.*y1;
t175 = t173.*y2;
t176 = t173.*y3;
t177 = -t173;
t178 = -t174;
t179 = -t175;
t180 = -t176;
t181 = t155+t178;
t182 = t156+t179;
t183 = t157+t180;
FF_med_uon = [y4;y5;y6;t60.*y1+t159.*y11;t60.*y2+t159.*y12;t60.*y3+t159.*y13;-U2.*U4.*t26.*t27.*t150.*t152;t48.*y11+t80.*y1+t29.*t46.*t155+U1.*t25.*t27.*t28.*t152.*t181;t48.*y12+t80.*y2+t29.*t46.*t156+U1.*t25.*t27.*t28.*t152.*t182;t48.*y13+t80.*y3+t29.*t46.*t157+U1.*t25.*t27.*t28.*t152.*t183;-y8;-y9;-y10;-U4.*t14.*t24.*t26.*t30.*t46.*t150];
if nargout > 1
    mt1 = [0.0,0.0,0.0,t60+t2.*t8.*t49.*y1.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t167.*y11,t2.*t8.*t49.*y2.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t167.*y12,t2.*t8.*t49.*y3.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t167.*y13,-U2.*U4.*t26.*t27.*t152.*t167+U2.*U4.*t26.*t27.*t150.*t153.*t169];
    mt2 = [t80+t160-t49.*y1.*y11.*3.0-t2.*t8.*t49.*y11.*3.0+t2.*t8.*t39.*t50.*y1.*5.0+U1.*t25.*t27.*t28.*t152.*(t154+t177+U4.*t14.*t24.*t26.*t47.*t81.*t148.*y1+U4.*t14.*t24.*t26.*t47.*t84.*t161.*y1+U4.*t2.*t8.*t14.*t24.*t26.*t60.*t84.*t148.*y1-U4.*t14.*t24.*t26.*t47.*t81.*t149.*t150.*t152.*y1-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t163.*y1-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t167.*y1+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t169.*y1+U4.*t2.*t8.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y1)-U1.*t25.*t27.*t28.*t153.*t169.*t181+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t81.*t148.*y1+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t161.*y1+U4.*t2.*t8.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y1];
    mt3 = [t49.*y2.*y11.*-3.0-t2.*t8.*t49.*y12.*3.0+t2.*t8.*t39.*t50.*y2.*5.0+U1.*t25.*t27.*t28.*t152.*(U4.*t14.*t24.*t26.*t47.*t81.*t148.*y2+U4.*t14.*t24.*t26.*t47.*t84.*t161.*y2+U4.*t2.*t8.*t14.*t24.*t26.*t60.*t84.*t148.*y2-U4.*t14.*t24.*t26.*t47.*t81.*t149.*t150.*t152.*y2-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t163.*y2-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t167.*y2+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t169.*y2+U4.*t2.*t8.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y2)-U1.*t25.*t27.*t28.*t153.*t169.*t182+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t81.*t148.*y2+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t161.*y2+U4.*t2.*t8.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y2];
    mt4 = [t49.*y3.*y11.*-3.0-t2.*t8.*t49.*y13.*3.0+t2.*t8.*t39.*t50.*y3.*5.0+U1.*t25.*t27.*t28.*t152.*(U4.*t14.*t24.*t26.*t47.*t81.*t148.*y3+U4.*t14.*t24.*t26.*t47.*t84.*t161.*y3+U4.*t2.*t8.*t14.*t24.*t26.*t60.*t84.*t148.*y3-U4.*t14.*t24.*t26.*t47.*t81.*t149.*t150.*t152.*y3-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t163.*y3-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t167.*y3+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t169.*y3+U4.*t2.*t8.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y3)-U1.*t25.*t27.*t28.*t153.*t169.*t183+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t81.*t148.*y3+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t161.*y3+U4.*t2.*t8.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y3,0.0,0.0,0.0,-U4.*t14.*t24.*t26.*t30.*t46.*t167,0.0,0.0,0.0,t3.*t9.*t49.*y1.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t168.*y11];
    mt5 = [t60+t3.*t9.*t49.*y2.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t168.*y12,t3.*t9.*t49.*y3.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t168.*y13,-U2.*U4.*t26.*t27.*t152.*t168+U2.*U4.*t26.*t27.*t150.*t153.*t171];
    mt6 = [t49.*y1.*y12.*-3.0-t3.*t9.*t49.*y11.*3.0+t3.*t9.*t39.*t50.*y1.*5.0+U1.*t25.*t27.*t28.*t152.*(U4.*t14.*t24.*t26.*t47.*t82.*t148.*y1+U4.*t14.*t24.*t26.*t47.*t84.*t162.*y1+U4.*t3.*t9.*t14.*t24.*t26.*t60.*t84.*t148.*y1-U4.*t14.*t24.*t26.*t47.*t82.*t149.*t150.*t152.*y1-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t165.*y1-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t168.*y1+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t171.*y1+U4.*t3.*t9.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y1)-U1.*t25.*t27.*t28.*t153.*t171.*t181+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t82.*t148.*y1+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t162.*y1+U4.*t3.*t9.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y1];
    mt7 = [t80+t160-t49.*y2.*y12.*3.0-t3.*t9.*t49.*y12.*3.0+t3.*t9.*t39.*t50.*y2.*5.0+U1.*t25.*t27.*t28.*t152.*(t154+t177+U4.*t14.*t24.*t26.*t47.*t82.*t148.*y2+U4.*t14.*t24.*t26.*t47.*t84.*t162.*y2+U4.*t3.*t9.*t14.*t24.*t26.*t60.*t84.*t148.*y2-U4.*t14.*t24.*t26.*t47.*t82.*t149.*t150.*t152.*y2-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t165.*y2-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t168.*y2+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t171.*y2+U4.*t3.*t9.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y2)-U1.*t25.*t27.*t28.*t153.*t171.*t182+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t82.*t148.*y2+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t162.*y2+U4.*t3.*t9.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y2];
    mt8 = [t49.*y3.*y12.*-3.0-t3.*t9.*t49.*y13.*3.0+t3.*t9.*t39.*t50.*y3.*5.0+U1.*t25.*t27.*t28.*t152.*(U4.*t14.*t24.*t26.*t47.*t82.*t148.*y3+U4.*t14.*t24.*t26.*t47.*t84.*t162.*y3+U4.*t3.*t9.*t14.*t24.*t26.*t60.*t84.*t148.*y3-U4.*t14.*t24.*t26.*t47.*t82.*t149.*t150.*t152.*y3-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t165.*y3-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t168.*y3+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t171.*y3+U4.*t3.*t9.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y3)-U1.*t25.*t27.*t28.*t153.*t171.*t183+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t82.*t148.*y3+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t162.*y3+U4.*t3.*t9.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y3,0.0,0.0,0.0,-U4.*t14.*t24.*t26.*t30.*t46.*t168,0.0,0.0,0.0,t4.*t10.*t49.*y1.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t170.*y11];
    mt9 = [t4.*t10.*t49.*y2.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t170.*y12,t60+t4.*t10.*t49.*y3.*3.0-U4.*t14.*t24.*t26.*t29.*t51.*t170.*y13,-U2.*U4.*t26.*t27.*t152.*t170+U2.*U4.*t26.*t27.*t150.*t153.*t172];
    mt10 = [t49.*y1.*y13.*-3.0-t4.*t10.*t49.*y11.*3.0+t4.*t10.*t39.*t50.*y1.*5.0+U1.*t25.*t27.*t28.*t152.*(U4.*t14.*t24.*t26.*t47.*t83.*t148.*y1+U4.*t14.*t24.*t26.*t47.*t84.*t164.*y1+U4.*t4.*t10.*t14.*t24.*t26.*t60.*t84.*t148.*y1-U4.*t14.*t24.*t26.*t47.*t83.*t149.*t150.*t152.*y1-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t166.*y1-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t170.*y1+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t172.*y1+U4.*t4.*t10.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y1)-U1.*t25.*t27.*t28.*t153.*t172.*t181+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t83.*t148.*y1+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t164.*y1+U4.*t4.*t10.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y1];
    mt11 = [t49.*y2.*y13.*-3.0-t4.*t10.*t49.*y12.*3.0+t4.*t10.*t39.*t50.*y2.*5.0+U1.*t25.*t27.*t28.*t152.*(U4.*t14.*t24.*t26.*t47.*t83.*t148.*y2+U4.*t14.*t24.*t26.*t47.*t84.*t164.*y2+U4.*t4.*t10.*t14.*t24.*t26.*t60.*t84.*t148.*y2-U4.*t14.*t24.*t26.*t47.*t83.*t149.*t150.*t152.*y2-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t166.*y2-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t170.*y2+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t172.*y2+U4.*t4.*t10.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y2)-U1.*t25.*t27.*t28.*t153.*t172.*t182+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t83.*t148.*y2+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t164.*y2+U4.*t4.*t10.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y2];
    mt12 = [t80+t160-t49.*y3.*y13.*3.0-t4.*t10.*t49.*y13.*3.0+t4.*t10.*t39.*t50.*y3.*5.0+U1.*t25.*t27.*t28.*t152.*(t154+t177+U4.*t14.*t24.*t26.*t47.*t83.*t148.*y3+U4.*t14.*t24.*t26.*t47.*t84.*t164.*y3+U4.*t4.*t10.*t14.*t24.*t26.*t60.*t84.*t148.*y3-U4.*t14.*t24.*t26.*t47.*t83.*t149.*t150.*t152.*y3-U4.*t14.*t24.*t26.*t47.*t84.*t150.*t152.*t166.*y3-U4.*t14.*t24.*t26.*t47.*t84.*t149.*t152.*t170.*y3+U4.*t14.*t24.*t26.*t47.*t84.*t149.*t150.*t153.*t172.*y3+U4.*t4.*t10.*t14.*t24.*t26.*t48.*t84.*t149.*t150.*t152.*y3)-U1.*t25.*t27.*t28.*t153.*t172.*t183+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t83.*t148.*y3+U4.*t14.*t24.*t26.*t29.*t46.*t47.*t84.*t164.*y3+U4.*t4.*t10.*t14.*t24.*t26.*t29.*t46.*t60.*t84.*t148.*y3,0.0,0.0,0.0,-U4.*t14.*t24.*t26.*t30.*t46.*t170,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0];
    mt13 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,U4.*t14.*t24.*t26.*t30.*t51.*t150.*y11,U4.*t14.*t24.*t26.*t30.*t51.*t150.*y12,U4.*t14.*t24.*t26.*t30.*t51.*t150.*y13,0.0,-t30.*t46.*t155,-t30.*t46.*t156,-t30.*t46.*t157,0.0,0.0,0.0,U4.*t14.*t24.*t26.*t29.^3.*t46.*t150.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,t159+U4.*t5.*t11.*t14.*t24.*t26.*t29.*t52.*t150.*y11,U4.*t5.*t11.*t14.*t24.*t26.*t29.*t52.*t150.*y12,U4.*t5.*t11.*t14.*t24.*t26.*t29.*t52.*t150.*y13,0.0,t48-t49.*y1.^2.*3.0+t5.*t11.*t29.*t51.*t155,t64+t5.*t11.*t29.*t51.*t156];
    mt14 = [t65+t5.*t11.*t29.*t51.*t157,0.0,0.0,0.0,-U4.*t5.*t11.*t14.*t24.*t26.*t30.*t51.*t150,0.0,0.0,0.0,U4.*t6.*t12.*t14.*t24.*t26.*t29.*t52.*t150.*y11,t159+U4.*t6.*t12.*t14.*t24.*t26.*t29.*t52.*t150.*y12,U4.*t6.*t12.*t14.*t24.*t26.*t29.*t52.*t150.*y13,0.0,t64+t6.*t12.*t29.*t51.*t155,t48-t49.*y2.^2.*3.0+t6.*t12.*t29.*t51.*t156,t66+t6.*t12.*t29.*t51.*t157,0.0,0.0,0.0,-U4.*t6.*t12.*t14.*t24.*t26.*t30.*t51.*t150,0.0,0.0,0.0,U4.*t7.*t13.*t14.*t24.*t26.*t29.*t52.*t150.*y11,U4.*t7.*t13.*t14.*t24.*t26.*t29.*t52.*t150.*y12,t159+U4.*t7.*t13.*t14.*t24.*t26.*t29.*t52.*t150.*y13,0.0,t65+t7.*t13.*t29.*t51.*t155,t66+t7.*t13.*t29.*t51.*t156,t48-t49.*y3.^2.*3.0+t7.*t13.*t29.*t51.*t157,0.0,0.0,0.0,-U4.*t7.*t13.*t14.*t24.*t26.*t30.*t51.*t150,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    mt15 = [U1.*t25.*t27.*t152.*t181,U1.*t25.*t27.*t152.*t182,U1.*t25.*t27.*t152.*t183,0.0,0.0,0.0,0.0];
    A_med_uon = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12,mt13,mt14,mt15],14,14);
end
end

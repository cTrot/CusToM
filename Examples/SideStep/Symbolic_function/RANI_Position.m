function out1 = RANI_Position(in1,in2,in3)
%RANI_POSITION
%    OUT1 = RANI_POSITION(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    25-Nov-2020 09:54:38

R1cut1_1 = in3(1);
R1cut1_2 = in3(4);
R1cut1_3 = in3(7);
R1cut2_1 = in3(2);
R1cut2_2 = in3(5);
R1cut2_3 = in3(8);
R1cut3_1 = in3(3);
R1cut3_2 = in3(6);
R1cut3_3 = in3(9);
p1cut1 = in2(1);
p1cut2 = in2(2);
p1cut3 = in2(3);
q1 = in1(1,:);
q23 = in1(23,:);
q24 = in1(24,:);
q25 = in1(25,:);
q26 = in1(26,:);
t2 = cos(q1);
t3 = cos(q23);
t4 = cos(q24);
t5 = cos(q25);
t6 = cos(q26);
t7 = sin(q1);
t8 = sin(q23);
t9 = sin(q24);
t10 = sin(q25);
t11 = sin(q26);
t12 = R1cut1_1.*t2;
t13 = R1cut1_2.*t2;
t14 = R1cut2_1.*t2;
t15 = R1cut2_2.*t2;
t16 = R1cut3_1.*t2;
t17 = R1cut3_2.*t2;
t18 = R1cut1_3.*t4;
t19 = R1cut2_3.*t4;
t20 = R1cut3_3.*t4;
t21 = R1cut1_1.*t7;
t22 = R1cut1_2.*t7;
t23 = R1cut2_1.*t7;
t24 = R1cut2_2.*t7;
t25 = R1cut3_1.*t7;
t26 = R1cut3_2.*t7;
t27 = R1cut1_3.*t9;
t28 = R1cut2_3.*t9;
t29 = R1cut3_3.*t9;
t30 = -t18;
t31 = -t19;
t32 = -t20;
t33 = -t21;
t34 = -t23;
t35 = -t25;
t36 = t12+t22;
t37 = t14+t24;
t38 = t16+t26;
t39 = t13+t33;
t40 = t15+t34;
t41 = t17+t35;
t42 = t3.*t36;
t43 = t3.*t37;
t44 = t3.*t38;
t45 = t8.*t36;
t46 = t8.*t37;
t47 = t8.*t38;
t48 = t3.*t39;
t49 = t3.*t40;
t50 = t3.*t41;
t51 = t8.*t39;
t52 = -t45;
t53 = t8.*t40;
t54 = -t46;
t55 = t8.*t41;
t56 = -t47;
t57 = t42+t51;
t58 = t43+t53;
t59 = t44+t55;
t60 = t48+t52;
t61 = t49+t54;
t62 = t50+t56;
t66 = -t4.*(t45-t48);
t67 = -t4.*(t46-t49);
t68 = -t4.*(t47-t50);
t69 = -t9.*(t45-t48);
t70 = -t9.*(t46-t49);
t71 = -t9.*(t47-t50);
t78 = -t10.*(t18+t9.*(t45-t48));
t79 = -t10.*(t19+t9.*(t46-t49));
t80 = -t10.*(t20+t9.*(t47-t50));
t63 = t5.*t57;
t64 = t5.*t58;
t65 = t5.*t59;
t72 = t27+t66;
t73 = t28+t67;
t74 = t29+t68;
t75 = t30+t69;
t76 = t31+t70;
t77 = t32+t71;
t81 = t63+t78;
t82 = t64+t79;
t83 = t65+t80;
out1 = [R1cut1_3.*(1.7e+1./2.25e+2)+p1cut1-t12.*7.555555555555556e-3-t13.*(1.7e+1./3.6e+2)+t21.*(1.7e+1./3.6e+2)-t22.*7.555555555555556e-3-t27.*4.297222222222222e-1-t63.*(1.7e+1./7.2e+2)-t5.*(t18+t9.*(t45-t48)).*(1.7e+1./6.0e+2)+t10.*(t18+t9.*(t45-t48)).*(1.7e+1./7.2e+2)-t10.*t57.*(1.7e+1./6.0e+2)-t6.*t72.*4.232244444444444e-1+t11.*t72.*(1.7e+1./9.0e+2)+t6.*t81.*(1.7e+1./9.0e+2)+t11.*t81.*4.232244444444444e-1+t4.*(t45-t48).*4.297222222222222e-1;R1cut2_3.*(1.7e+1./2.25e+2)+p1cut2-t14.*7.555555555555556e-3-t15.*(1.7e+1./3.6e+2)+t23.*(1.7e+1./3.6e+2)-t24.*7.555555555555556e-3-t28.*4.297222222222222e-1-t64.*(1.7e+1./7.2e+2)-t5.*(t19+t9.*(t46-t49)).*(1.7e+1./6.0e+2)+t10.*(t19+t9.*(t46-t49)).*(1.7e+1./7.2e+2)-t10.*t58.*(1.7e+1./6.0e+2)-t6.*t73.*4.232244444444444e-1+t11.*t73.*(1.7e+1./9.0e+2)+t6.*t82.*(1.7e+1./9.0e+2)+t11.*t82.*4.232244444444444e-1+t4.*(t46-t49).*4.297222222222222e-1;R1cut3_3.*(1.7e+1./2.25e+2)+p1cut3-t16.*7.555555555555556e-3-t17.*(1.7e+1./3.6e+2)+t25.*(1.7e+1./3.6e+2)-t26.*7.555555555555556e-3-t29.*4.297222222222222e-1-t65.*(1.7e+1./7.2e+2)-t5.*(t20+t9.*(t47-t50)).*(1.7e+1./6.0e+2)+t10.*(t20+t9.*(t47-t50)).*(1.7e+1./7.2e+2)-t10.*t59.*(1.7e+1./6.0e+2)-t6.*t74.*4.232244444444444e-1+t11.*t74.*(1.7e+1./9.0e+2)+t6.*t83.*(1.7e+1./9.0e+2)+t11.*t83.*4.232244444444444e-1+t4.*(t47-t50).*4.297222222222222e-1];
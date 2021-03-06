function dx1dF1 = fse_dx1dF1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1)
%FSE_DX1DF1
%    DX1DF1 = FSE_DX1DF1(F1,F2,A1,A2,A3,A4,A7,ALPHA,BETA,X1)

%    This function was generated by the Symbolic Math Toolbox version 5.6.
%    23-Aug-2013 13:40:26

t3 = F2.*a1;
t2 = F1-t3;
t14 = alpha.*t2;
t4 = beta-t14;
t5 = exp(t4);
t6 = t5+1.0;
t7 = 1.0./t6;
t8 = a7+1.0;
t9 = 1.0./a3;
t21 = t9.*x1;
t10 = -t21+1.0;
t11 = 1.0./t8;
t12 = t10.^t11;
t13 = a2-1.0;
t15 = a3-x1;
t33 = alpha.*t15;
t16 = -beta-t33;
t17 = exp(t16);
t18 = t17+1.0;
t19 = 1.0./t18;
t20 = t13.*t2.*t7;
t22 = a3.*t12.*t8;
t23 = t20+t22;
t30 = alpha.*t23;
t24 = beta-t30;
t25 = exp(t24);
t26 = t25+1.0;
t27 = 1.0./t26;
t28 = t11.*t13.*t2.*t7.*t9;
t29 = t12+t28;
t31 = t29.^t8;
t32 = 1.0./t6.^2;
t34 = a3.*t27.*t31;
t35 = -a3+t34+x1;
t36 = t2.*t9.*x1;
t37 = t19.*t35.*t7;
t38 = a3+t37-x1;
t43 = alpha.*t38;
t39 = beta-t43;
t40 = exp(t39);
t41 = t40+1.0;
t42 = 1.0./t41;
t44 = t11.*t13.*t7.*t9;
t45 = alpha.*t11.*t13.*t2.*t32.*t5.*t9;
t46 = t44+t45;
t47 = t29.^a7;
t48 = a3.*t27.*t46.*t47.*t8;
t49 = 1.0./t26.^2;
t50 = t13.*t7;
t51 = alpha.*t13.*t2.*t32.*t5;
t52 = t50+t51;
t53 = a3.*alpha.*t25.*t31.*t49.*t52;
t54 = t48+t53;
t55 = t19.*t54.*t7;
t56 = alpha.*t15.*t32.*t5;
t57 = alpha.*t19.*t32.*t35.*t5;
t58 = t55+t56+t57;
t59 = t15.*t7;
t60 = t37+t59;
t61 = a4.*t9.*x1.*(1.0./2.0);
t62 = t36+t61;
t63 = a3.*t7;
t71 = t42.*t60;
t64 = t63-t71;
t70 = t62.*t7;
t72 = a4.*t64.*t9.*(1.0./2.0);
t65 = t36-t70-t72+x1;
t76 = alpha.*t65;
t66 = -beta-t76;
t67 = exp(t66);
t68 = t67+1.0;
t69 = 1.0./t68;
t73 = t55+t57;
t74 = 1.0./t41.^2;
t75 = alpha.*t40.*t60.*t73.*t74;
t77 = t7-1.0;
t78 = t42.*t58;
t79 = a3.*alpha.*t32.*t5;
dx1dF1 = -t75+t79-t42.*t58-t69.*t77.*(t21+a4.*t9.*(t75+t78-a3.*alpha.*t32.*t5).*(1.0./2.0)-t7.*t9.*x1-alpha.*t32.*t5.*t62)-alpha.*t32.*t5.*t65.*t69-alpha.*t65.*t67.*1.0./t68.^2.*t77.*(t21+a4.*t9.*(t75+t78-t79).*(1.0./2.0)-t7.*t9.*x1-alpha.*t32.*t5.*t62);

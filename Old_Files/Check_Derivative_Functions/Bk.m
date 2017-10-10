function [B] = Bk(alpha,beta,F,a,x)
a1 = a(1); a2 = a(2); a3 = a(3); a4 = a(4); a5 = a(5); a6 = a(6); a7 = a(7);

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5);

F1 = F(1); F2 = F(2);

alpha = alpha.alpha1;
beta = beta.beta1;

B(1,1) = fse_dx1da1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1); 
B(1,2) = fse_dx1da2(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(1,3) = fse_dx1da3(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(1,4) = fse_dx1da4(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(1,5) = 0;
B(1,6) = 0;
B(1,7) = fse_dx1da7(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);

B(2,1) = fse_dx2da1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(2,2) = fse_dx2da2(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(2,3) = fse_dx2da3(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(2,4) = fse_dx2da4(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(2,5) = 0;
B(2,6) = -x2;
B(2,7) = fse_dx2da7(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);

B(3,1) = fse_dx3da1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(3,2) = fse_dx3da2(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(3,3) = fse_dx3da3(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(3,4) = fse_dx3da4(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
B(3,5) = -x3;
B(3,6) = 0;
B(3,7) = fse_dx3da7(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);

B(4,1) = 0;
B(4,2) = 0;
B(4,3) = 0;
B(4,4) = 0; 
B(4,5) = x3-x4;
B(4,6) = 0;
B(4,7) = 0;

B(5,1) = 0;
B(5,2) = 0;
B(5,3) = 0;
B(5,4) = 0;
B(5,5) = x4-x5;
B(5,6) = 0;
B(5,7) = 0;
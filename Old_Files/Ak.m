function [A] = Ak(alpha,beta,F,a,x)
a1 = a(1); a2 = a(2); a3 = a(3); a4 = a(4); a5 = a(5); a6 = a(6); a7 = a(7);

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5);

F1 = F(1); F2 = F(2);

% alpha1 = alpha.alpha1; alpha2 = alpha.alpha2; alpha3 = alpha.alpha3; alpha5 = alpha.alpha5; alpha6 = alpha.alpha6; alpha7 = alpha.alpha7; alpha9 = alpha.alpha9; alpha10 = alpha.alpha10; alpha11 = alpha.alpha11; alpha12 = alpha.alpha12;
% beta1 = beta.beta1; beta3 = beta.beta3; beta9 = beta.beta9; beta10 = beta.beta10; beta11 = beta.beta11; beta12 = beta.beta12;

alpha = alpha.alpha1;
beta = beta.beta1;

A(1,1) = fse_dx1dx1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1); 
A(1,2) = 0; 
A(1,3) = 0; 
A(1,4) = 0; 
A(1,5) = 0;

A(2,1) = fse_dx2dx1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
A(2,2) = -a6+1.0;
A(2,3) = 0;
A(2,4) = 0;
A(2,5) = 0;

A(3,1) = fse_dx3dx1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1);
A(3,2) = 0;
A(3,3) = -a5+1.0;
A(3,4) = 0;
A(3,5) = 0;

A(4,1) = 0;
A(4,2) = 0;
A(4,3) = a5;
A(4,4) = -a5+1.0;
A(4,5) = 0;

A(5,1) = 0;
A(5,2) = 0;
A(5,3) = 0;
A(5,4) = a5;
A(5,5) = -a5+1.0;
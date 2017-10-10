clear; 
clc;
%Adding constraints and computing first derivatives

syms x1 x2 x3 x4 x5 a1 a2 a3 a4 a5 a6 a7 F1 F2 alpha beta lambda1 lambda2;

%Water Balance Computations
%Soil Precipitation
cond1 = (F1-a1*F2); %P - aET   
approx_f1 = 1+exp(-alpha*(cond1)+beta);
M1 = ((F1-a1*F2)*(1-a2))/approx_f1;

%Lagrangian model
LM1 = M1 - (lambda1*a1 + lambda2*a2);

%Derivatives
%Original Model
dM1da1 = diff(M1,a1);
matlabFunction(dM1da1, 'file', 'fse_dM1da1', 'vars', symvar(dM1da1));

dM1da2 = diff(M1,a2);
matlabFunction(dM1da2, 'file', 'fse_dM1da2', 'vars', symvar(dM1da2));

%Lagrangian Model
dLM1da1 = diff(LM1,a1);
% matlabFunction(dLM1da1, 'file', 'fse_dLM1da1', 'vars', symvar(dLM1da1));

dLM1da2 = diff(LM1,a2);
% matlabFunction(dLM1da2, 'file', 'fse_dLM1da2', 'vars', symvar(dLM1da2));

%% Synthetic Modeling
ref_PP  = dlmread('/Users/humbertovergara/Documents/Trabajo/Parameter_and_QPEs_Error_Characterization_for_DA/TS_Data/MPE_tarboro_basin_avg_time_series.csv', ',', 1,1);
PP = ref_PP(:,2);
PP_period = datenum('01-Jan-2002 01:00:00'):1/24:datenum('31-Dec-2009 23:00:00');

PET = dlmread('~/Documents/MS_Work/Thesis_Work/HyMOD/PET_mmhr_Tarboro_2002-2009.txt', ' ', 1,0);
PET_period = datenum('01-Jan-2002 00:00:00'):1/24:datenum('31-Dec-2009 23:00:00');

%"Real Values"
Cet = 0.5;
Imp = 0.3;

cd1 = (PP(1:100)-Cet*PET(1:100)); %P - aET   
app1 = 1+exp(-1000*(cd1)+100);
obs = ((PP(1:100)-Cet.*PET(1:100)).*(1-Imp))./app1;

%Incorrect first estimates
Cet = 0.9;
Imp = 0.9;

cd1 = (PP(1:100)-Cet*PET(1:100)); %P - aET   
app1 = 1+exp(-1000*(cd1)+100);
model = ((PP(1:100)-Cet.*PET(1:100)).*(1-Imp))./app1;

%% FSM
% syms U H;
%Number of Variables
n = 0; %States
p = 2; %Parameters

%Jacobian of h
% Dh = zeros(1,n,nsteps);
Dh = zeros(1,100);

%Set V(0) and U(0)
% V(:,:,1) = eye(n);
% U(:,:,1) = zeros(n,p);

V(1) = 0;
U = zeros(p,100);

for k = 1:100
        %Jacobian of h
        Dh(1,k) = 1;

        %Hessian 
%         H(k,:) = [Dh(:,:,k)*V(:,:,k), Dh(:,:,k)*U(:,:,k)];
        H(k,:) = Dh(k)*U(:,k);

        %Compute A(k)
        %compute CCREST's forward sensitivity equations w.r.t. initial conditions
%         A(:,:,k) = Ak(alpha,beta,F(:,k),a,x(:,k));
        A(k) = 0;
        
        %Compute B(k)
        %compute CCREST's forward sensitivity equations w.r.t. model parameters
%         B(:,:,k) = Bk(alpha,beta,F(:,k),a,x(:,k));
%         B(1,k) = fse_dLM1da1(PP(k),PET(k),Cet,Imp,1000,100);
%         B(2,k) = fse_dLM1da2(PP(k),PET(k),Cet,1000,100);
        
        B(1,k) = fse_dM1da1(PP(k),PET(k),Cet,Imp,1000,100);
        B(2,k) = fse_dM1da2(PP(k),PET(k),Cet,1000,100);

        %Compute V(k+1)
%         V(:,:,k+1) = A(:,:,k)*V(:,:,k);
%         V(k+1) = A(k)*V(k);

        %Compute U(k+1)
%         U(:,:,k+1) = A(:,:,k)*U(:,:,k)+B(:,:,k);
        U(:,k+1) = B(:,k);
        
        %Compute the forecast error
        ef(k) = obs(k) - model(k);
end
ef(isnan(ef)==1) = 0;

% if (n+p < length(obs))
%     %Over-determined
%     condnum = cond(H'*H);
%     rankcoef = rank(H'*H);
    correction = inv(H'*H)*H'*ef';
% else
%     %Under-determined
% %     condnum = cond(H'*H);
% %     rankcoef = rank(H'*H);
%     correction = (H'/(H*H'))*ef';
% end

%Minimize cost function for a1, a2, lambda1, and lambda2
cost_func = 0.5*((ef' - H*correction)'*1*(ef' - H*correction));


Cet = Cet+correction(1);
Imp = Imp+correction(2);

cd1 = (PP(1:100)-Cet*PET(1:100)); %P - aET   
app1 = 1+exp(-1000*(cd1)+100);
new_model = ((PP(1:100)-Cet.*PET(1:100)).*(1-Imp))./app1;

plot(obs, 'DisplayName', 'Obs');
hold all;
plot(model, 'DisplayName', 'Incorrect');
plot(new_model, 'DisplayName', 'Correct');

% %Infiltration
% cond2 = (a3*(1+a7)*(1-(x1/a3))^(1/(1+a7))-M1);
% approx_f2 = 1 + exp(-alpha*cond2+beta);
% cond3 = (a3-x1);
% approx_f3 = 1 + exp(-alpha*cond3-beta);
% M2 = (a3 - x1 - (a3*((1-(x1/a3))^(1/(1+a7))-(M1/(a3*(1+a7))))^(1+a7))/approx_f2)*(1/(approx_f3*approx_f1));
% L_M2 = M2 - lambda*x1;
% 
% %EXPRESSIONS FOR DERIVATIVES (Sensitivities)
% 
% %Derivatives w.r.t Forcing?
% dLM2dx1 = diff(L_M2,x1);
% matlabFunction(dLM2dx1, 'file', 'fse_dLM2dx1', 'vars', symvar(dLM2dx1));

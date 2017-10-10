function [correction, V, U, A, H, Dh, condnum,rankcoef] = fsmccrest(rundata,pred,alpha,beta)

%"Unpack" rundata
Pin = rundata.PP;
PET = rundata.PET;
% States_0 = rundata.States_0;
Parameters = rundata.Parameters;

%Simulation Settings
tperiod = rundata.tperiod; %time period for simulation
basin_area = rundata.barea; %basin area in km^2
nsteps = length(tperiod);

%Subset data for simulation
st_idx = find(rundata.PET_period == tperiod(1));
en_idx = find(rundata.PET_period == tperiod(end));
PET = PET(st_idx:en_idx);

st_idx = find(rundata.PP_period == tperiod(1));
en_idx = find(rundata.PP_period == tperiod(end));
Pin = Pin(st_idx:en_idx);

st_idx = find(rundata.obsQ_period == rundata.tperiod(1));
en_idx = find(rundata.obsQ_period == rundata.tperiod(end));
obsQ = rundata.obsQ(st_idx:en_idx);

%Forcing
F(1,:) = Pin; %Precipitation (mm)
F(2,:) = PET; %Potential Evapotranspiration (mm)

%Model Parameters
a(1) = Parameters(1); %ET Factor
a(2) = Parameters(2)/100; %Impervious Area Percentage
a(3) = Parameters(3); %Soil's Maximum Water Capacity (mm)
a(4) = Parameters(4); %Ksat in (mm)
a(5) = Parameters(5); %Overland Reservoir Depletion Rate
a(6) = Parameters(6); %Interflow Reservoir Depletion Rate
a(7) = Parameters(7); %Exponential in VIC model

%Model States
x = zeros(5,length(F));
% x(1,:) = a(3).*(pred.SM(1:end-1)/100); %Soil Moisture(mm)
x(1,:) = pred.SM(1:end-1); %Soil Moisture(mm)
x(2,:) = pred.ISU(1:end-1); %Interflow Reservoir (mm)
x(3:end,:) = pred.ISO(:,1:end-1); %Overland Reservoirs (mm)

if (isreal(x) == 0)
	error('Imaginary Prediction');
end

%Number of Variables
n = 5; %States
p = 7; %Parameters

%Jacobian of h
Dh = zeros(1,n,nsteps);

%Set V(0) and U(0)
V(:,:,1) = eye(n);
U(:,:,1) = zeros(n,p);

%Initialize some variables
% H = zeros(nsteps,n+p);
% A = zeros(n,n,nsteps);
% B = zeros(n,p,nsteps);

% k = 1;
% cval = round(pred.NonLinear(:,1));
% for ki = 1:nsteps
for k = 1:nsteps
    %check for jump
%     jump = abs(round(pred.NonLinear(:,ki)) - cval);
%     if (max(jump) > 0)
        %Jacobian of h
        Dh(:,:,k) = jacob_ic(a);

        %Hessian 
        H(k,:) = [Dh(:,:,k)*V(:,:,k), Dh(:,:,k)*U(:,:,k)];

        %Compute A(k)
        %compute CCREST's forward sensitivity equations w.r.t. initial conditions
        A(:,:,k) = Ak(alpha,beta,F(:,k),a,x(:,k));
        
        %Compute B(k)
        %compute CCREST's forward sensitivity equations w.r.t. model parameters
        B(:,:,k) = Bk(alpha,beta,F(:,k),a,x(:,k));
        
        %Compute V(k+1)
        V(:,:,k+1) = A(:,:,k)*V(:,:,k);

        %Compute U(k+1)
        U(:,:,k+1) = A(:,:,k)*U(:,:,k)+B(:,:,k);
        
        %Compute the forecast error
        ef(k) = (obsQ(k) - pred.Q(k)).*(3.6/basin_area);
%         temp_Q(k) = pred.Q(ki);
%         k = k + 1;
%     end

        if (max(isnan(A(:))) > 0)
%             test = Bk(alpha,beta,F(:,k),a,x(:,k));
            test = Ak(alpha,beta,F(:,k),a,x(:,k));
            error('NaNs');
        end
end
ef(isnan(ef)==1) = 0;

if (n+p < length(obsQ))
    %Over-determined
    condnum = cond(H'*H);
    rankcoef = rank(H'*H);
    correction = pinv(H'*H)*H'*ef';
else
    %Under-determined
%     condnum = cond(H'*H);
%     rankcoef = rank(H'*H);
    correction = (H'/(H*H'))*ef';
end

% correction = []; condnum = []; rankcoef = [];

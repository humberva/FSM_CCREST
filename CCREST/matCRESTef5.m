function [cSM, aET, cERI, cERO, infiltration] = matCRESTef5(stepHours, precipIn, petIn, parameters, states)
                                 
precip = precipIn .* stepHours; % precipIn is mm/hr, precip is mm
pet = petIn .* stepHours; % petIn in mm/hr, pet is mm 

[nrows, ncols] = size(precip);

infiltration = zeros(nrows, ncols);
precipSoil = zeros(nrows, ncols);
precipImperv = zeros(nrows, ncols);
interflowExcess = zeros(nrows, ncols);

Wmaxm = zeros(nrows, ncols);
A = zeros(nrows, ncols);
R = zeros(nrows, ncols);
cERI = zeros(nrows, ncols);
cERO = zeros(nrows, ncols);
Wo = zeros(nrows, ncols);
temX = zeros(nrows, ncols);
ExcessET = zeros(nrows, ncols);

adjPET = pet * parameters.PKE;

%We have more water coming in than leaving via ET.
precipSoil(precip > adjPET) = (precip(precip > adjPET) - adjPET(precip > adjPET)) .* (1 - parameters.PIM(precip > adjPET)); %This is the precip that makes it to the soil
precipImperv(precip > adjPET) = precip(precip > adjPET) - adjPET(precip > adjPET) - precipSoil(precip > adjPET); %Portion of precip on impervious area

interflowExcess(precip > adjPET) = states.SM(precip > adjPET) - parameters.PWM(precip > adjPET);

interflowExcess(precip > adjPET & interflowExcess < 0) = 0.0;

states.SM(precip > adjPET & states.SM > parameters.PWM) = parameters.PWM(precip > adjPET & states.SM > parameters.PWM);

Wmaxm(precip > adjPET & states.SM < parameters.PWM) = parameters.PWM(precip > adjPET & states.SM < parameters.PWM) .* (1 + parameters.PB(precip > adjPET & states.SM < parameters.PWM));
A(precip > adjPET & states.SM < parameters.PWM) = Wmaxm(precip > adjPET & states.SM < parameters.PWM) .* (1 - (1 - states.SM(precip > adjPET & states.SM < parameters.PWM)  ./ parameters.PWM(precip > adjPET & states.SM < parameters.PWM)) .^ ( 1 ./ (1 + parameters.PB(precip > adjPET & states.SM < parameters.PWM))));
      
R(precip > adjPET & states.SM < parameters.PWM & precipSoil + A >= Wmaxm) = precipSoil(precip > adjPET & states.SM < parameters.PWM & precipSoil + A >= Wmaxm) - (parameters.PWM(precip > adjPET & states.SM < parameters.PWM & precipSoil + A >= Wmaxm) - states.SM(precip > adjPET & states.SM < parameters.PWM & precipSoil + A >= Wmaxm)); % Leftovers after filling SM

R(R < 0) = 0;

Wo(precip > adjPET & states.SM < parameters.PWM & precipSoil + A >= Wmaxm) = parameters.PWM(precip > adjPET & states.SM < parameters.PWM & precipSoil + A >= Wmaxm);

infiltration(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) = parameters.PWM(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) .* ((1 - A(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) ./ Wmaxm(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm)) .^ (1 + parameters.PB(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm)) - (1 - (A(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) + precipSoil(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm)) ./ Wmaxm(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm)).^(1 + parameters.PB(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm)));

infiltration(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm & infiltration > precipSoil) = precipSoil(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm & infiltration > precipSoil);

R(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) = precipSoil(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) - infiltration(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm);

R(R < 0) = 0;

Wo(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) = states.SM(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm) + infiltration(precip > adjPET & states.SM < parameters.PWM & precipSoil + A < Wmaxm);
               
R(precip > adjPET & states.SM >= parameters.PWM) = precipSoil(precip > adjPET & states.SM >= parameters.PWM);
Wo(precip > adjPET & states.SM >= parameters.PWM) = parameters.PWM(precip > adjPET & states.SM >= parameters.PWM);
infiltration(precip > adjPET & states.SM >= parameters.PWM) = 0.0; %Added by HV so we get the infiltration value as output

% Now R is excess water, split it between overland & interflow

temX(precip > adjPET) = ((states.SM(precip > adjPET) + Wo(precip > adjPET)) ./ (parameters.PWM(precip > adjPET) .* 2)) .* (parameters.PFC(precip > adjPET) .* stepHours); % Calculate how much water can infiltrate


cERI(precip > adjPET & R <= temX) = R(precip > adjPET & R <= temX);
cERI(precip > adjPET & R > temX) = temX(precip > adjPET & R > temX);

cERO(precip > adjPET) = R(precip > adjPET) - cERI(precip > adjPET) + precipImperv(precip > adjPET);

aET(precip > adjPET) = adjPET(precip > adjPET);

cERI(precip > adjPET) = cERI(precip > adjPET) + interflowExcess(precip > adjPET); % Extra interflow that got routed.

% All the incoming precip goes straight to ET

infiltration(precip <= adjPET) = 0.0; %Added by HV so we get the infiltration value as output
cERO(precip <= adjPET) = 0.0;
interflowExcess(precip <= adjPET) = states.SM(precip <= adjPET) - parameters.PWM(precip <= adjPET);

interflowExcess(precip <= adjPET & interflowExcess < 0) = 0.0;

cERI(precip <= adjPET) = interflowExcess(precip <= adjPET);

states.SM(precip <= adjPET & states.SM > parameters.PWM) = parameters.PWM(precip <= adjPET & states.SM > parameters.PWM);

ExcessET(precip <= adjPET) = (adjPET(precip <= adjPET) - precip(precip <= adjPET)) .* states.SM(precip <= adjPET) ./ parameters.PWM(precip <= adjPET);
Wo(precip <= adjPET & ExcessET < states.SM) = states.SM(precip <= adjPET & ExcessET < states.SM) - ExcessET(precip <= adjPET & ExcessET < states.SM); % We can evaporate away ExcessET too.
Wo(precip <= adjPET & ExcessET >= states.SM) = 0.0; % We don't have enough to evaporate ExcessET.
ExcessET(precip <= adjPET & ExcessET >= states.SM) = states.SM(precip <= adjPET & ExcessET >= states.SM);
aET(precip <= adjPET) = ExcessET(precip <= adjPET) + precip(precip <= adjPET);


%states.SM = Wo;
cSM = Wo; %New Soil moisture

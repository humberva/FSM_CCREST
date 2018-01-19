function [cSM, aET, cERI, cERO, infiltration] = CRESTef5(stepHours, precipIn, petIn, parameters, states)
                                 
precip = precipIn * stepHours; % precipIn is mm/hr, precip is mm
pet = petIn * stepHours; % petIn in mm/hr, pet is mm
infiltration = 0;

%R = 0.0; 
%Wo = 0.0;

adjPET = pet * parameters.PKE;
%temX = 0.0;

%We have more water coming in than leaving via ET.
if (precip > adjPET)
    precipSoil = (precip - adjPET) * (1 - parameters.PIM); %This is the precip that makes it to the soil
    precipImperv = precip - adjPET - precipSoil; %Portion of precip on impervious area

    interflowExcess = states.SM - parameters.PWM;
    if (interflowExcess < 0)
      interflowExcess = 0;
    end

    if (states.SM > parameters.PWM)
      states.SM = parameters.PWM;
    end

    if (states.SM < parameters.PWM)
        Wmaxm = parameters.PWM * (1 + parameters.PB);
        A = Wmaxm * (1 - (1 - states.SM  / parameters.PWM) ^ ( 1 / (1 + parameters.PB)));
      
        if (precipSoil + A >= Wmaxm)
            %infiltration = (parameters.PWM - states.SM), not explicit in EF5
            R = precipSoil - (parameters.PWM - states.SM); % Leftovers after filling SM

            if (R < 0)
                fprintf('R to %f, [%f, %f, %f, %f, %f, %f]\n', R, parameters.PWM, parameters.PB, states.SM, A, Wmaxm, precipSoil);
                R = 0;
            end

            Wo = parameters.PWM;
        else
            infiltration = parameters.PWM .* ((1 - A / Wmaxm) .^ (1 + parameters.PB) - (1 - (A + precipSoil) ./ Wmaxm).^(1 + parameters.PB));
%             infiltration = parameters.PWM - states.SM - parameters.PWM * ((1 - (A + precipSoil) / Wmaxm)^(1 + parameters.PB));
            
            
            if (infiltration > precipSoil)
                infiltration = precipSoil;
            else
                if (infiltration < 0)
                    fprintf('Infiltration went to %f, [%f, %f, %f, %f, %f, %f]\n', infiltration, parameters.PWM, parameters.PB, states.SM, A, Wmaxm, precipSoil);
                end
            end

            R = precipSoil - infiltration;

            if (R < 0)
                fprintf('R (infil) to %f, [%f, %f, %f, %f, %f, %f]\n', R, parameters.PWM, parameters.PB, states.SM, A, Wmaxm, precipSoil);
                R = 0;
            end
            Wo = states.SM + infiltration;
        end
        
    else
        R = precipSoil;
        Wo = parameters.PWM;
        infiltration = 0; %Added by HV so we get the infiltration value as output
    end

    % Now R is excess water, split it between overland & interflow

    temX = ((states.SM + Wo) / (parameters.PWM * 2)) * (parameters.PFC * stepHours); % Calculate how much water can infiltrate

    if (R <= temX)
      cERI = R;
    else
      cERI = temX;
    end
    
    cERO = R - cERI + precipImperv;

    aET = adjPET;

    cERI = cERI + interflowExcess; % Extra interflow that got routed.
else               % All the incoming precip goes straight to ET
    infiltration = 0; %Added by HV so we get the infiltration value as output
    cERO = 0;
    interflowExcess = states.SM - parameters.PWM;
    
    if (interflowExcess < 0)
      interflowExcess = 0;
    end
    
    cERI = interflowExcess;

    if (states.SM > parameters.PWM)
      states.SM = parameters.PWM;
    end

    ExcessET = (adjPET - precip) * states.SM / parameters.PWM;
    if (ExcessET < states.SM)
      Wo = states.SM - ExcessET; % We can evaporate away ExcessET too.
    else
      Wo = 0; % We don't have enough to evaporate ExcessET.
      ExcessET = states.SM;
    end
    aET = ExcessET + precip;
end

%states.SM = Wo;
cSM = Wo; %New Soil moisture

% Add Overland Excess Water to fastFlow
%  *fastFlow += (cERO / (stepHours * 3600.0f));
% Add Interflow Excess Water to slowFlow
%  *slowFlow += (cERI / (stepHours * 3600.0f));
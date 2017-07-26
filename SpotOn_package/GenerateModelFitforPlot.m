function [model_PDF, model_CDF] = GenerateModelFitforPlot(model_params, JumpProb, JumpProbCDF, NumberOfStates)
% GenerateModelPDFforPlot 
%   This function is for plotting the model output as both PDF histogram
%   and CDF

global LocError dT HistVecJumps dZ HistVecJumpsCDF ModelFit FitLocError Z_corr_a Z_corr_b JumpsPerdT UseWeights
Temp_ModelFit = ModelFit;

%%%%%%%%%%%%%%%%%%%%%%%% GENERATE MODEL-FIT PDF %%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Because ModelFit is a global variable, temporarily force-change it to
%   PDF-fitting (ModelFit=1) and then change it back to whatever it was
%   afterwards. 
ModelFit = 1; % force-change the global variable to PDF

% calculate the model PDF using the input model params:
if NumberOfStates == 2 && FitLocError == 0 % 2-state model, fixed Loc Error
    y = Model_2State(model_params, JumpProb);

elseif NumberOfStates == 2 && FitLocError == 1 % 2-state model, Loc Error from fitting
    y = Model_2State_fitLocError(model_params, JumpProb);
    
elseif NumberOfStates == 3 && FitLocError == 0 % 3-state model, fixed Loc Error
    y = Model_3State(model_params, JumpProb);
    
elseif NumberOfStates == 3 && FitLocError == 1 % 3-state model, Loc Error from fitting
    y = Model_3State_fitLocError(model_params, JumpProb);    
end

% Make model-output normalized PDF for plotting
model_PDF = zeros(size(y,1), size(y,2));
%Normalize y as a PDF
for i=1:size(y,1)
    model_PDF(i,:) = y(i,:)./sum(y(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% GENERATE MODEL-FIT CDF %%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Because ModelFit is a global variable, temporarily force-change it to
%   CDF-fitting (ModelFit=2) and then change it back to whatever it was
%   afterwards. 
ModelFit = 2; % force-change the global variable to CDF

% calculate the model PDF using the input model params:
if NumberOfStates == 2 && FitLocError == 0 % 2-state model, fixed Loc Error
    model_CDF = Model_2State(model_params, JumpProbCDF);

elseif NumberOfStates == 2 && FitLocError == 1 % 2-state model, Loc Error from fitting
    model_CDF = Model_2State_fitLocError(model_params, JumpProbCDF);
    
elseif NumberOfStates == 3 && FitLocError == 0 % 3-state model, fixed Loc Error
    model_CDF = Model_3State(model_params, JumpProbCDF);
    
elseif NumberOfStates == 3 && FitLocError == 1 % 3-state model, Loc Error from fitting
    model_CDF = Model_3State_fitLocError(model_params, JumpProbCDF);    
end

if UseWeights == 1
    %Normalize CDF to get rid of weighting
    for i=1:length(JumpsPerdT)
        model_CDF(i,:) = model_CDF(i,:)./JumpsPerdT(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Force change the global variable back to whatever it was at the end of
% this function to avoid causing problems:
ModelFit = Temp_ModelFit;
end


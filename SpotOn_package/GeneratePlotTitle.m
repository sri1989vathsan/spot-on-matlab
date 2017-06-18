function [ PlotTitle ] = GeneratePlotTitle(SampleName, NumberOfStates, model_params, Min3Traj, CellLocs, CellJumps, CellFrames, TrajNumb)
%GeneratePlotTitle makes a title for a plot
%   makes titles for both single-cell and merged datasets

global FitLocError ModelFit

% main title
main_name = SampleName;

% model-specific names:
if NumberOfStates == 2 && FitLocError == 0 % 2-state model, fixed Loc Error
    model_specific_name = {['Dfree = ', num2str(round(model_params(1)*1000)/1000), ' um2/s; Dbound = ', num2str(round(model_params(2)*1000)/1000), ...
    ' um2/s; FracBound = ', num2str(round(model_params(3)*1000)/1000)];};
    
elseif NumberOfStates == 2 && FitLocError == 1 % 2-state model, Loc Error from fitting
    model_specific_name = {['Dfree = ', num2str(round(model_params(1)*1000)/1000), ' um2/s; Dbound = ', num2str(round(model_params(2)*1000)/1000), ' um2/s;']; ...
        ['FracBound = ', num2str(round(model_params(3)*1000)/1000), '; fitted LocError = ', num2str(round(model_params(4)*1000)/1000), ' um'];};
    
elseif NumberOfStates == 3 && FitLocError == 0 % 3-state model, fixed Loc Error
    model_specific_name = {['Dfree1 = ', num2str(round(model_params(1)*1000)/1000), ' um2/s; FracFree1 = ', num2str(round(model_params(5)*1000)/1000)];...
        ['Dfree2 = ', num2str(round(model_params(2)*1000)/1000), ' um2/s; FracFree2 = ', num2str(round((1-model_params(4)-model_params(5))*1000)/1000)];...
        ['Dbound = ', num2str(round(model_params(3)*1000)/1000), ' um2/s; FracBound = ', num2str(round(model_params(4)*1000)/1000)];  };
    
elseif NumberOfStates == 3 && FitLocError == 1 % 3-state model, Loc Error from fitting
    model_specific_name = {['Dfree1 = ', num2str(round(model_params(1)*1000)/1000), ' um2/s; FracFree1 = ', num2str(round(model_params(5)*1000)/1000)];...
        ['Dfree2 = ', num2str(round(model_params(2)*1000)/1000), ' um2/s; FracFree2 = ', num2str(round((1-model_params(4)-model_params(5))*1000)/1000)];...
        ['Dbound = ', num2str(round(model_params(3)*1000)/1000), ' um2/s; FracBound = ', num2str(round(model_params(4)*1000)/1000), '; fitted LocError = ', num2str(round(model_params(6)*1000)/1000), ' um'];  };
end

% output PDF of CDF fitting:
if ModelFit == 1
    fit_type = 'PDF';
elseif ModelFit == 2
    fit_type = 'CDF';
end

% overview of number of localizations etc. 
localizations_overview_name = {['Total trajs: ', num2str(TrajNumb), '; trajs >=3: ', num2str(Min3Traj), '; FitType: ', fit_type];...
 ['Locs = ', num2str(CellLocs), '; Locs/Frame = ', num2str(round(CellLocs/CellFrames*1000)/1000), ...
'; jumps: ', num2str(CellJumps), ';']};

% now merge all the titles:
PlotTitle = [main_name; model_specific_name; localizations_overview_name;];
end


function [Output_struct] = SpotOn_core(Params)
%SPOTON_CORE core SpotOn function 
%   The function takes the user-input and performs data processing and
%   function fitting

% define neccesary global parameters:
global LocError dT HistVecJumps dZ HistVecJumpsCDF ModelFit FitLocError Z_corr_a Z_corr_b
colour = jet; close;
Output_struct = struct([]); % for saving outputs

%%%%%%%%%%%%%%%%% Unpack the input structure Params: %%%%%%%%%%%%%%%%%%%%%%
TimeGap = Params.TimeGap; 
dZ = Params.dZ; 
dT = Params.TimeGap/1000; % convert to seconds
GapsAllowed = Params.GapsAllowed; 
TimePoints = Params.TimePoints; 
BinWidth = Params.BinWidth; 
UseAllTraj = Params.UseAllTraj;
JumpsToConsider = Params.JumpsToConsider; 
MaxJumpPlotPDF = Params.MaxJumpPlotPDF; 
MaxJumpPlotCDF = Params.MaxJumpPlotCDF; 
MaxJump = Params.MaxJump;
SavePlot = Params.SavePlot; 
ModelFit = Params.ModelFit;
DoSingleCellFit = Params.DoSingleCellFit; 
FitIterations = Params.FitIterations; 
FitLocErrorRange = Params.FitLocErrorRange; 
LocError = Params.LocError; 
FitLocError = Params.FitLocError;
NumberOfStates = Params.NumberOfStates;
D_Free_2State = Params.D_Free_2State; 
D_Bound_2State = Params.D_Bound_2State; 
D_Free1_3State = Params.D_Free1_3State; 
D_Free2_3State = Params.D_Free2_3State; 
D_Bound_3State = Params.D_Bound_3State;
curr_dir = Params.curr_dir; 
data_struct = Params.data_struct; 
SampleName = Params.SampleName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% If saving plots, make sure that the folder exists %%%%%%%%%%%%%
save_path = [curr_dir, filesep, 'SavedPlots', filesep];
if ~exist(save_path, 'dir')
    % the directory is not there, so make it
    mkdir(save_path)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Define some additional variables %%%%%%%%%%%%%%%%%%%%%%%
HistVecJumps = 0:BinWidth:MaxJump; % histogram/PDF displacement bins in micrometers
HistVecJumpsCDF = 0:0.001:MaxJump; % CDF displacement bins in micrometers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% PARAMETERS FOR Z-CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%
% find the correct parameters for z-correction; these have been
% pre-computed through Monte Carlo simulations and will be inferred by
% matching dT and dZ to their closest values from a library.
[Z_corr_a, Z_corr_b] = MatchZ_corr_coeff(dT, dZ, GapsAllowed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PART 1: Analyze data from each single cell %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DoSingleCellFit == 1 
    %   Analyze each single cell individually
    %   1. Load data for each single cell
    %   2. Make histrogram of jump lengths
    %   3. Fit model to each single cell
    %   4. Plot the model fit and best fit parameters
    
    % initialize various counters
    last_integer = 0;
    PlotIndex = 0;
    CellNumb = 1;

    % loop over all replicates and all cells in each replicate
    for RepIter=1:length(data_struct)
        workspaces = data_struct(RepIter).workspaces;
        curr_path = data_struct(RepIter).path;
        for WorkIter=1:length(workspaces)
        
            %Can only plot 8 cells at a time, so create a new figure if
            %neccesary:
            PlotIndex = PlotIndex + 1;
            if ceil(CellNumb/8) > last_integer
                %Plot everything on a big figure:
                figure('position',[100 100 1600 1200]); %[x y width height]
                last_integer = ceil(CellNumb/8);
                PlotIndex = 1;
            end
            
            
            %%%%% STEP 1: LOAD IN A CELL AND COMPILE HISTOGRAMS
            tic; disp(['Analysing cell ', num2str(WorkIter), ' of ', num2str(length(workspaces)), ' in replicate ', num2str(RepIter), ' of ', num2str(length(data_struct))]);
            % use function "compile_histograms_single_cell.m" to compile histograms
            disp('loading in trajectories and compiling histograms...');
            %full_path = [curr_path, char(workspaces(WorkIter)), '.mat'] % full_path of workspace to be loaded
            [JumpProb, JumpProbCDF, Min3Traj, CellLocs, CellJumps, CellJumps_used, CellFrames, TrajNumb] = compile_histograms_single_cell([curr_path, char(workspaces(WorkIter))], UseAllTraj, GapsAllowed, TimePoints, JumpsToConsider);
            
            
            %%%%% STEP 2: PERFORM MODEL-FITTING OF THE CURRENT CELL
            disp('performing model fitting of displacement histograms...');
            [model_params, ~, ~] = ModelFitting_main(JumpProb, JumpProbCDF, NumberOfStates, FitLocErrorRange, FitIterations, D_Free_2State, D_Bound_2State, D_Free1_3State, D_Free2_3State, D_Bound_3State);
            toc;

            %%%%% STEP 3: PLOT THE BEST FIT FOR THE CURRENT CELL
            % get a normalized PDF/CDF histogram for plotting using current
            % best-fit parameters:
            disp('proceeding to plotting the output of the model fitting...');
            [model_PDF, model_CDF] = GenerateModelFitforPlot(model_params, JumpProb, JumpProbCDF, NumberOfStates);

            % Generate a plot title with all the relevant info
            PlotTitle = GeneratePlotTitle(workspaces{RepIter}, NumberOfStates, model_params, Min3Traj, CellLocs, CellJumps, CellJumps_used, CellFrames, TrajNumb);
            
            % Do actual plotting
            subplot(2,4,PlotIndex);
            histogram_spacer = 0.055;
            hold on;
            for i=size(JumpProb,1):-1:1
                new_level = (i-1)*histogram_spacer;%*y_max;
                colour_element = colour(round(i/size(JumpProb,1)*size(colour,1)),:);
                plot(HistVecJumps, new_level*ones(1,length(HistVecJumps)), 'k-', 'LineWidth', 1); 
                for j=2:size(JumpProb,2)
                    x1 = HistVecJumps(1,j-1); x2 = HistVecJumps(1,j);
                    y1 = new_level; y2 = JumpProb(i,j-1)+new_level;
                    patch([x1 x1 x2 x2], [y1 y2 y2 y1],colour_element);
                end
                plot(HistVecJumps, model_PDF(i,:)+new_level, 'k-', 'LineWidth', 2);
                text(0.6*MaxJumpPlotPDF, new_level+0.5*histogram_spacer, ['\Delta\tau: ', num2str(TimeGap*i), ' ms'], 'HorizontalAlignment','left', 'FontSize',9, 'FontName', 'Helvetica');
            end
            axis([0 MaxJumpPlotPDF 0 1.05*(max(JumpProb(end,:))+(size(JumpProb,1)-1)*histogram_spacer)]);
            title(PlotTitle, 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
            set(gca,'YColor','w')
            ylabel('Probability', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
            xlabel('jump length \mu m', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
            set(gca,'YTickLabel',''); set(gca, 'YTick', []);
            hold off;

            % save all the relevant info from fitting this cell:
            if CellNumb == 1
                single_cell_model_params = model_params;
                single_cell_model_PDF = model_PDF;
                single_cell_model_CDF = model_CDF;
                single_cell_data_PDF = JumpProb;
                single_cell_data_CDF = JumpProbCDF;
            else
                single_cell_model_params = vertcat(single_cell_model_params, model_params);
                single_cell_model_PDF = vertcat(model_PDF, single_cell_model_PDF);
                single_cell_model_CDF = vertcat(model_CDF, single_cell_model_CDF);
                single_cell_data_PDF = vertcat(single_cell_data_PDF, JumpProb);
                single_cell_data_CDF = vertcat(single_cell_data_CDF, JumpProbCDF);
            end
            
            % Increment CellNumb:
            CellNumb = CellNumb+1;
        end
    end
    % save the single-cell plots:
    if SavePlot == 1
        print([save_path, SampleName, '_SingleCell_FitType',num2str(ModelFit), 'FitLocError', num2str(FitLocError),'.pdf'], '-dpdf');
        print([save_path, SampleName, '_SingleCell_FitType',num2str(ModelFit), 'FitLocError', num2str(FitLocError),'.eps'], '-depsc');
    end
    Output_struct(1).single_cell_model_params = single_cell_model_params;
    Output_struct(1).single_cell_model_PDF = single_cell_model_PDF;
    Output_struct(1).single_cell_model_CDF = single_cell_model_CDF;
    Output_struct(1).single_cell_data_PDF = single_cell_data_PDF;
    Output_struct(1).single_cell_data_CDF = single_cell_data_CDF;
    clear JumpProb JumpProbCDF model_params 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PART 2: Analyze data from all cells %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OVERVIEW:
    %   1. Load data and merged data for all cells
    %   2. Make histrograms (PDF and CDF) of displacements
    %   3. Fit model to the merged dara
    %   4. Plot the model fit for both CDF and PDF and best fit parameters

% load in merged data and compile histograms:
disp('===========================================================');
disp('merging data from all cells and replicates...');
disp('loading in trajectories and compiling histograms...'); tic;
%%%%% STEP 1: LOAD IN MERGED DATA AND COMPILE HISTOGRAMS
[JumpProb, JumpProbCDF, Min3Traj, TotalLocs, TotalFrames, TotalJumps, TotalJumps_used, TrajNumb, DyeSurvivalProb, DyeHistVec, DyeMean] = compile_histograms_many_cells( data_struct, UseAllTraj, GapsAllowed, TimePoints, JumpsToConsider );
toc;
%%%%% STEP 2: PERFORM MODEL-FITTING OF THE MERGED DATA
disp('performing model fitting of displacement histograms...'); tic;
[model_params, residuals, list_of_model_parameters] = ModelFitting_main(JumpProb, JumpProbCDF, NumberOfStates, FitLocErrorRange, FitIterations, D_Free_2State, D_Bound_2State, D_Free1_3State, D_Free2_3State, D_Bound_3State);
toc;
%%%%% STEP 3: PLOT THE BEST FIT FOR THE MERGED DATA
% get a normalized PDF/CDF histogram for plotting using current
% best-fit parameters:
disp('proceeding to plotting the output of the model fitting...');
[model_PDF, model_CDF] = GenerateModelFitforPlot(model_params, JumpProb, JumpProbCDF, NumberOfStates);
% Generate a plot title with all the relevant info
PlotTitle = GeneratePlotTitle(SampleName, NumberOfStates, model_params, Min3Traj, TotalLocs, TotalJumps, TotalJumps_used, TotalFrames, TrajNumb);
       
% PLOT THE SURVIVAL PROBABILITY OF THE FLUOROPHORE
figure('position',[800 100 300 275]); %[x y width height]
hold on;
plot(DyeHistVec, DyeSurvivalProb, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
axis([1 51 0.001 1.01]);
title(['1-CDF of trajectory lengths; mean = ', num2str(DyeMean), ' frames'], 'FontSize',9, 'FontName', 'Helvetica');
ylabel('1-CDF', 'FontSize',9, 'FontName', 'Helvetica');
xlabel('number of frames', 'FontSize',9, 'FontName', 'Helvetica');
set(gca,'yscale','log');
hold off;


% PLOT THE RESIDUALS FROM THE FITTING
% Plot residuals for the relevant fit: 
%   so PDF residuals for PDF-fitting
%   or CDF residuals for CDF-fitting
figure('position',[300 300 1200 800]); %[x y width height]
% find max_y for plot
max_y = max([0.1 max(max(abs(residuals)))]); min_y = -max_y;
for i=1:min([12 size(residuals,1)])
    colour_element = colour(round(i/size(residuals,1)*size(colour,1)),:);
    subplot(3,4,i);
    hold on;
    if ModelFit == 1
        plot(HistVecJumps, residuals(i,:), '-', 'Color', colour_element, 'LineWidth', 2);
        max_x = MaxJumpPlotPDF;
    elseif ModelFit == 2
        plot(HistVecJumpsCDF, residuals(i,:), '-', 'Color', colour_element, 'LineWidth', 2);
        max_x = MaxJumpPlotCDF;
    end
    plot([0 max_x], [0 0], 'k--', 'LineWidth', 1);
    
    axis([0 max_x min_y max_y]);
    if ModelFit == 1
        title({SampleName; ['PDF residuals for \Delta\tau: ', num2str(TimeGap*i), ' ms']}, 'FontSize',8, 'FontName', 'Helvetica', 'Color', 'k');
    elseif ModelFit == 2
        title({SampleName; ['CDF residuals for \Delta\tau: ', num2str(TimeGap*i), ' ms']}, 'FontSize',8, 'FontName', 'Helvetica', 'Color', 'k');
    end
    ylabel('residuals', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
    xlabel('displacements (\mu m)', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
    hold off;
end
if SavePlot == 1
    print([save_path, SampleName, '_residuals_FitType',num2str(ModelFit), '.pdf'], '-dpdf');
    print([save_path, SampleName, '_residuals_FitType',num2str(ModelFit), '.eps'], '-depsc');
end

%PLOT CDFs of DISPLACEMENTS AND OF FIT
figure('position',[100 100 1200 800]); %[x y width height]
for i=1:min([12 size(JumpProbCDF,1)])
    colour_element = colour(round(i/size(JumpProbCDF,1)*size(colour,1)),:);
    subplot(3,4,i);
    hold on;
    plot(HistVecJumpsCDF, JumpProbCDF(i,:), '-', 'LineWidth', 2, 'Color', colour_element);
    plot(HistVecJumpsCDF, model_CDF(i,:), 'k-', 'LineWidth', 1);
    
    axis([0 MaxJumpPlotCDF 0 1.05]);
    title({SampleName; ['CDF for \Delta\tau: ', num2str(TimeGap*i), ' ms']}, 'FontSize',8, 'FontName', 'Helvetica', 'Color', 'k');
    ylabel('displacement CDF', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
    xlabel('displacements (\mu m)', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
    legend('raw data', 'Model fit', 'Location', 'SouthEast');
    legend boxoff
    hold off;
end
if SavePlot == 1
    print([save_path, SampleName, '_mergedCDFs_FitType',num2str(ModelFit), '.pdf'], '-dpdf');
    print([save_path, SampleName, '_mergedCDFs_FitType',num2str(ModelFit), '.eps'], '-depsc');
end

%PLOT THE HISTOGRAM OF TRANSLOCATIONS
figure('position',[200 200 300 400]); %[x y width height]
histogram_spacer = 0.055;
hold on;
for i=size(JumpProb,1):-1:1
    new_level = (i-1)*histogram_spacer;
    colour_element = colour(round(i/size(JumpProb,1)*size(colour,1)),:);
    plot(HistVecJumps, new_level*ones(1,length(HistVecJumps)), 'k-', 'LineWidth', 1); 
    for j=2:size(JumpProb,2)
        x1 = HistVecJumps(1,j-1); x2 = HistVecJumps(1,j);
        y1 = new_level; y2 = JumpProb(i,j-1)+new_level;
        patch([x1 x1 x2 x2], [y1 y2 y2 y1],colour_element);
    end
    plot(HistVecJumps, model_PDF(i,:)+new_level, 'k-', 'LineWidth', 2);
    text(0.6*MaxJumpPlotPDF, new_level+0.5*histogram_spacer, ['\Delta\tau: ', num2str(TimeGap*i), ' ms'], 'HorizontalAlignment','left', 'FontSize',9, 'FontName', 'Helvetica');
end
axis([0 MaxJumpPlotPDF 0 1.05*(max(JumpProb(end,:))+(size(JumpProb,1)-1)*histogram_spacer)]);
title(PlotTitle, 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
set(gca,'YColor','w')
ylabel('Probability', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
xlabel('jump length \mu m', 'FontSize',9, 'FontName', 'Helvetica', 'Color', 'k');
set(gca,'YTickLabel',''); set(gca, 'YTick', []);
hold off;
if SavePlot == 1
    print([save_path, SampleName, '_mergedPDFs_FitType',num2str(ModelFit), 'FitLocError', num2str(FitLocError), '.pdf'], '-dpdf');
    print([save_path, SampleName, '_mergedPDFs_FitType',num2str(ModelFit), 'FitLocError', num2str(FitLocError),'.eps'], '-depsc');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%%%%%%%%%%%%%%%%%%%%%%%% SAVE THE DESIRED OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%
% save key model outputs for all merged cells:
Output_struct(1).merged_model_params = model_params;
Output_struct(1).merged_model_PDF = model_PDF;
Output_struct(1).merged_model_CDF = model_CDF;
Output_struct(1).merged_data_PDF = JumpProb;
Output_struct(1).merged_data_CDF = JumpProbCDF;
Output_struct(1).merged_model_residuals = residuals;
Output_struct(1).list_of_model_parameters = list_of_model_parameters;
Output_struct(1).r_bins_PDF = HistVecJumps;
Output_struct(1).r_bins_CDF = HistVecJumpsCDF;

end




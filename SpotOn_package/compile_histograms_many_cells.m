function [JumpProb, JumpProbCDF, Min3Traj, TotalLocs, TotalFrames, TotalJumps, TotalJumps_used, TrajNumb, DyeSurvivalProb, DyeHistVec, DyeMean] = compile_histograms_many_cells( data_struct, UseAllTraj, GapsAllowed, TimePoints, JumpsToConsider )
%compile_histograms_many_cells Compile histograms for many cells
%   the input will be a structure containing the neccesary info to load in
%   trackedPar from many single cells
global HistVecJumps HistVecJumpsCDF

%%%%% initialize key data-storage parameters
AllData = []; %for storing data
TotalFrames = 0; %for counting how many frames
TotalLocs = 0; %for counting the total number of localizations
%POOLED DATA: input is a strucutred array
%Load data
for iter=1:length(data_struct)
    currPath = data_struct(iter).path;
    currWorkspaces = data_struct(iter).workspaces;
    currInclude = data_struct(iter).Include;
    % Check that Workspaces and Include lengths match
    if exist('currWorkspaces')
        if length(currWorkspaces) ~= length(currInclude)
            error(['Workspace vector and Include vector do not have the same length - fix before running; Workspace length = ', num2str(length(workspaces)), '; Include length = ', num2str(length(Include)), ';'])
        end
    end
    for i=1:length(currWorkspaces)
        if currInclude(i) == i
            full_path = [currPath, char(currWorkspaces(i))];
            % check to see whether the user added ".mat" to the end of the workspace
            % MAT file:
            if full_path(end-3:end) == '.mat'
                load(full_path, 'trackedPar'); % Load in trackedPar from the relevant cell
            else
                load([full_path, '.mat'], 'trackedPar'); % Load in trackedPar from the relevant cell
            end
            AllData = [AllData trackedPar]; 
            for n=1:length(trackedPar)
                TotalLocs = TotalLocs + length(trackedPar(1,n).Frame);
            end
            %Find total frames using a slightly ad-hoc way: find the last frame with
            %a localization and round it. This is not an elegant solution, but it
            %works for most reasonable particle densities:
            LastIdx = length(trackedPar);
            TempLastFrame = max(trackedPar(1,LastIdx).Frame);
            TotalFrames = TotalFrames + 100*round(TempLastFrame/100);
        end
    end
end

Min3Traj = 0; %for counting number of min3 trajectories;
TotalJumps = 0; %for counting the total number of jumps
TotalJumps_used = 0; %for counting the total number of jumps actually used
TrajLengthHist = zeros(1,length(AllData));
TrajNumb = length(AllData); % total number of trajectories

%Calculate a histogram of translocation lengths vs. frames
TransFrames = TimePoints+GapsAllowed*(TimePoints-1); TransLengths = struct; 
for i=1:TransFrames
    TransLengths(1,i).Step = []; %each iteration is a different number of timepoints
end
if UseAllTraj == 1 %Use all of the trajectory
    for i=1:length(AllData)
        % save length of the trajectory
        TrajLengthHist(1,i) = max(AllData(i).Frame) - min(AllData(i).Frame) + 1;
        CurrTrajLength = size(AllData(i).xy,1);
        %save lengths
        if CurrTrajLength >= 3
            Min3Traj = Min3Traj + 1;
        end
        %Now loop through the trajectory. Keep in mind that there are missing
        %timepoints in the trajectory, so some gaps may be for multiple
        %timepoints.
        %Figure out what the max jump to consider is:
        HowManyFrames = min(TimePoints-1, CurrTrajLength);
        if CurrTrajLength > 1
            TotalJumps = TotalJumps + CurrTrajLength - 1; %for counting all the jumps
            for n=1:HowManyFrames
                for k=1:CurrTrajLength-n
                    %Find the current XY coordinate and frames between
                    %timepoints
                    CurrXY_points = vertcat(AllData(i).xy(k,:), AllData(i).xy(k+n,:));
                    CurrFrameJump = AllData(i).Frame(k+n) - AllData(i).Frame(k);
                    %Calculate the distance between the pair of points
                    TransLengths(1,CurrFrameJump).Step = horzcat(TransLengths(1,CurrFrameJump).Step, pdist(CurrXY_points));
                end
            end
        end    
    end
    TotalJumps_used = TotalJumps; % all jumps were used, so these are the same
elseif UseAllTraj == 0 %Use only the first JumpsToConsider displacements
    for i=1:length(AllData)
        % save length of the trajectory
        TrajLengthHist(1,i) = max(AllData(i).Frame) - min(AllData(i).Frame) + 1;
        CurrTrajLength = size(AllData(i).xy,1);
        if CurrTrajLength >= 3
            Min3Traj = Min3Traj + 1;
        end
        %Loop through the trajectory. If it is a short trajectory, you need to
        %make sure that you do not overshoot. So first figure out how many
        %jumps you can consider.
        %Figure out what the max jump to consider is:
        HowManyFrames = min([TimePoints-1, CurrTrajLength]);
        if CurrTrajLength > 1
            TotalJumps = TotalJumps + CurrTrajLength - 1; %for counting all the jumps
            TotalJumps_used = TotalJumps_used + min([CurrTrajLength-1 JumpsToConsider]); %for counting all the jumps actually used
            for n=1:HowManyFrames
                FrameToStop = min([CurrTrajLength, n+JumpsToConsider]);
                for k=1:(FrameToStop-n)
                    %Find the current XY coordinate and frames between
                    %timepoints
                    CurrXY_points = vertcat(AllData(i).xy(k,:), AllData(i).xy(k+n,:));
                    CurrFrameJump = AllData(i).Frame(k+n) - AllData(i).Frame(k);
                    %Calculate the distance between the pair of points
                    TransLengths(1,CurrFrameJump).Step = horzcat(TransLengths(1,CurrFrameJump).Step, pdist(CurrXY_points));
                end
            end
        end  
    end
end    
% CALCULATE THE SURVIVAL PROBABILITY OF THE FLUOROPHORE
DyeHistVec = 1:1:max(TrajLengthHist);
TrajLengthProb = histc(TrajLengthHist, DyeHistVec)./length(TrajLengthHist);
TrajLengthCDF = zeros(1,length(TrajLengthProb));
for i=2:length(TrajLengthProb)
    TrajLengthCDF(1,i) = sum(TrajLengthProb(1:i));
end
DyeSurvivalProb = 1-TrajLengthCDF;
DyeMean = mean(TrajLengthHist);

% CALCULATE THE PDF HISTOGRAMS
JumpProb = zeros(TimePoints-1, length(HistVecJumps));
JumpProbFine = zeros(TimePoints-1, length(HistVecJumpsCDF));
for i=1:size(JumpProb,1)
    JumpProb(i,:) = histc(TransLengths(1,i).Step, HistVecJumps)/length(TransLengths(1,i).Step);
    JumpProbFine(i,:) = histc(TransLengths(1,i).Step, HistVecJumpsCDF)/length(TransLengths(1,i).Step);
end  

% CALCULATE THE CDF HISTOGRAMS:
JumpProbCDF = zeros(TimePoints-1, length(HistVecJumpsCDF));
for i=1:size(JumpProbCDF,1)
    for j=2:size(JumpProbCDF,2)
        JumpProbCDF(i,j) = sum(JumpProbFine(i,1:j));
    end
end  

end


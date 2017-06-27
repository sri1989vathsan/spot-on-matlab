function [JumpProb, JumpProbCDF, Min3Traj, CellLocs, CellJumps, CellJumps_used, CellFrames, TrajNumb] = compile_histograms_single_cell( full_path, UseAllTraj, GapsAllowed, TimePoints, JumpsToConsider )
%COMPILE_HISTOGRAMS_SINGLE_CELL Compiles histograms from a single cell
%   Use the function for compiling histograms and other relevant info from
%   a single cell
global HistVecJumps HistVecJumpsCDF

% check to see whether the user added ".mat" to the end of the workspace
% MAT file:
if full_path(end-3:end) == '.mat'
    load(full_path, 'trackedPar'); % Load in trackedPar from the relevant cell
else
    load([full_path, '.mat'], 'trackedPar'); % Load in trackedPar from the relevant cell
end

%Find total frames using a slight ad-hoc way: find the last frame with
%a localization and round it. This is not an elegant solution, but it
%works for most reasonable particle densities:
CellLocs = 0; % for counting the total number of localizations
LastIdx = length(trackedPar);
TempLastFrame = max(trackedPar(1,LastIdx).Frame);
CellFrames = 100*round(TempLastFrame/100);
for n=1:length(trackedPar)
    CellLocs = CellLocs + length(trackedPar(1,n).Frame);
end

% total number of trajectories
TrajNumb = length(trackedPar);

%Compile histograms for each jump length
Min3Traj = 0; %for counting number of min3 trajectories;
CellJumps = 0; %for counting the total number of jumps
CellJumps_used = 0; %for counting the total number of jumps actually used
TransFrames = TimePoints+GapsAllowed*(TimePoints-1); TransLengths = struct; 
for i=1:TransFrames
    TransLengths(1,i).Step = []; %each iteration is a different number of timepoints
end

if UseAllTraj == 1 %Use all of the trajectory
    for i=1:length(trackedPar)
        CurrTrajLength = size(trackedPar(i).xy,1);
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
            CellJumps = CellJumps + CurrTrajLength - 1; %for counting all the jumps
            for n=1:HowManyFrames
                for k=1:CurrTrajLength-n
                    %Find the current XY coordinate and frames between
                    %timepoints
                    CurrXY_points = vertcat(trackedPar(i).xy(k,:), trackedPar(i).xy(k+n,:));
                    CurrFrameJump = trackedPar(i).Frame(k+n) - trackedPar(i).Frame(k);
                    %Calculate the distance between the pair of points
                    TransLengths(1,CurrFrameJump).Step = horzcat(TransLengths(1,CurrFrameJump).Step, pdist(CurrXY_points));
                end
            end
        end    
    end
    CellJumps_used = CellJumps; % all jumps were used, so these are the same
elseif UseAllTraj == 0 % Use only the first JumpsToConsider displacements
    for i=1:length(trackedPar)
        CurrTrajLength = size(trackedPar(i).xy,1);
        if CurrTrajLength >= 3
            Min3Traj = Min3Traj + 1;
        end
        %Loop through the trajectory. If it is a short trajectory, you need to
        %make sure that you do not overshoot. So first figure out how many
        %jumps you can consider.
        %Figure out what the max jump to consider is:
        HowManyFrames = min([TimePoints-1, CurrTrajLength]);
        if CurrTrajLength > 1
            CellJumps = CellJumps + CurrTrajLength - 1; %for counting all the jumps
            CellJumps_used = CellJumps_used + min([CurrTrajLength-1 JumpsToConsider]); %for counting all the jumps actually used
            for n=1:HowManyFrames
                FrameToStop = min([CurrTrajLength, n+JumpsToConsider]);
                for k=1:(FrameToStop-n)
                    %Find the current XY coordinate and frames between
                    %timepoints
                    CurrXY_points = vertcat(trackedPar(i).xy(k,:), trackedPar(i).xy(k+n,:));
                    CurrFrameJump = trackedPar(i).Frame(k+n) - trackedPar(i).Frame(k);
                    %Calculate the distance between the pair of points
                    TransLengths(1,CurrFrameJump).Step = horzcat(TransLengths(1,CurrFrameJump).Step, pdist(CurrXY_points));
                end
            end
        end  
    end
end 

%CALCULATE THE PDF HISTOGRAMS
JumpProb = zeros(TimePoints-1, length(HistVecJumps));
JumpProbFine = zeros(TimePoints-1, length(HistVecJumpsCDF));
for i=1:size(JumpProb,1)
    JumpProb(i,:) = histc(TransLengths(1,i).Step, HistVecJumps)/length(TransLengths(1,i).Step);
    JumpProbFine(i,:) = histc(TransLengths(1,i).Step, HistVecJumpsCDF)/length(TransLengths(1,i).Step);
end  

%CALCULATE THE CDF HISTOGRAMS:
JumpProbCDF = zeros(TimePoints-1, length(HistVecJumpsCDF));
for i=1:size(JumpProbCDF,1)
    for j=2:size(JumpProbCDF,2)
        JumpProbCDF(i,j) = sum(JumpProbFine(i,1:j));
    end
end  
        


end


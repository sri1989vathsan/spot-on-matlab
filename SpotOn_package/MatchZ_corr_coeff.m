function [ Z_corr_a, Z_corr_b ] = MatchZ_corr_coeff( dT, dZ, GapsAllowed )
%MATCHZ_CORR_PARAMS 
%   Match the closest a,b coefficients for the Z-correction.
%   Use kd-tree to efficiently search the 2-dimensional space for the
%   coefficients for the Z-correction.
%   These coefficients were generated from Monte Carlo simulations: please
%   see the SpotOn documentation for full details. 
%   Code for generating kd-tree: Mdl = KDTreeSearcher(matrix_dTdZ);

% load in the database of a,b parameters:
if GapsAllowed == 0
    load('MonteCarloParams_0_gaps.mat', 'Mdl', 'a', 'b', 'matrix_dTdZ');
elseif GapsAllowed == 1
    load('MonteCarloParams_1_gap.mat', 'Mdl', 'a', 'b', 'matrix_dTdZ');
elseif GapsAllowed == 2
    load('MonteCarloParams_2_gaps.mat', 'Mdl', 'a', 'b', 'matrix_dTdZ');
elseif GapsAllowed > 2
    error(['user supplied GapsAllowed = ', num2str(GapsAllowed), '; but the Matlab version of Spot-On does not currently support more than 2 gaps.']);
end

% find the nearest match (Euclidean distance)
IdxNN = knnsearch(Mdl, [dT, dZ], 'K', 1);
Z_corr_a = a(IdxNN);
Z_corr_b = b(IdxNN);
matched_dT = matrix_dTdZ(IdxNN,1);
matched_dZ = matrix_dTdZ(IdxNN,2);

% output to screen which parameters were found
disp('============ finding Z-correction coeffcients ===========');
disp(['User supplied time gap: ', num2str(dT*1000), ' ms; matched time gap: ', num2str(matched_dT*1000), ' ms;']);
disp(['User supplied dZ observation slice: ', num2str(dZ*1000), ' nm; matched dZ: ', num2str(matched_dZ*1000), ' nm;']);
disp('Using the following coefficients:');
disp(['Z_corr_a = ', num2str(Z_corr_a), '; Z_corr_b = ', num2str(Z_corr_b)]);
disp('===== found coefficients: return to data processing =====');
end


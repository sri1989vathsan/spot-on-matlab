function Binned_y = Model_3State_fitLocError( parameter_guess, JumpProb )
%Model_3State_fitLocError Summary of this function goes here
%   Detailed explanation goes here

global dT HistVecJumps dZ HistVecJumpsCDF ModelFit Z_corr_a Z_corr_b JumpsPerdT UseWeights

% define key parameters
r = HistVecJumpsCDF;
y = zeros(size(JumpProb,1), length(r));
Binned_y_PDF = zeros(size(JumpProb,1), size(JumpProb,2));


%Assign a value to each parameter
D_FREE1 = parameter_guess(1);
D_FREE2 = parameter_guess(2);
D_BOUND = parameter_guess(3);
F_BOUND = parameter_guess(4);
F_FREE1 = parameter_guess(5);
FIT_LocError = parameter_guess(6);

%Assume ABSORBING BOUNDARIES
Z_corr1 = zeros(1,size(JumpProb,1)); Z_corr2 = zeros(1,size(JumpProb,1));

for iterator=1:size(JumpProb,1)
    %Calculate the jump length distribution of the parameters for each time-jump
    curr_dT = iterator*dT;
    
    %Calculate the axial Z-correction
    %First calculate the corrected DeltaZ:
    DeltaZ_use1 = dZ + Z_corr_a  * sqrt(D_FREE1) + Z_corr_b;
    DeltaZ_use2 = dZ + Z_corr_a  * sqrt(D_FREE2) + Z_corr_b;
    
    % use half DeltaZ
    HalfDeltaZ_use1 = DeltaZ_use1/2;
    HalfDeltaZ_use2 = DeltaZ_use2/2;
    
    %Compute the integral for each of the free componentd
    Z_corr1(1,iterator) =  1/DeltaZ_use1 * integral(@(z)C_AbsorBoundAUTO(z,curr_dT, D_FREE1, HalfDeltaZ_use1),-HalfDeltaZ_use1,HalfDeltaZ_use1);
    Z_corr2(1,iterator) =  1/DeltaZ_use2 * integral(@(z)C_AbsorBoundAUTO(z,curr_dT, D_FREE2, HalfDeltaZ_use2),-HalfDeltaZ_use2,HalfDeltaZ_use2);
    
    %update the function output
    y(iterator,:) = Z_corr1(1,iterator).*F_FREE1.*(r./(2*(D_FREE1*curr_dT+FIT_LocError^2))).*exp(-r.^2./(4*(D_FREE1*curr_dT+FIT_LocError^2))) ...
                    + Z_corr2(1,iterator).*(1-F_BOUND-F_FREE1).*(r./(2*(D_FREE2*curr_dT+FIT_LocError^2))).*exp(-r.^2./(4*(D_FREE2*curr_dT+FIT_LocError^2))) ...
                    + F_BOUND.*(r./(2*(D_BOUND*curr_dT+FIT_LocError^2))).*exp(-r.^2./(4*(D_BOUND*curr_dT+FIT_LocError^2))) ;
end

% Make sure the model output (PDF or CDF) fits the model input:
if ModelFit == 1
    % Now bin the output y so that it matches the JumpProb variable: 
    for i=1:size(JumpProb,1)
        for j=1:size(JumpProb,2)
            if j == size(JumpProb,2)
                Binned_y_PDF(i,j) = mean(y(i,maxIndex:end));
            else
                [~, minIndex] = min(abs(r-HistVecJumps(j)));
                [~, maxIndex] = min(abs(r-HistVecJumps(j+1)));
                Binned_y_PDF(i,j) = mean(y(i,minIndex:maxIndex-1));
            end
        end
    end

    %Normalize:
    for i=1:size(JumpProb,1)
        Binned_y_PDF(i,:) = Binned_y_PDF(i,:)./sum(Binned_y_PDF(i,:));
    end

    %You want to fit to a histogram
    %So no need to calculate the CDF
    Binned_y = Binned_y_PDF;
    
    if UseWeights == 1
        % Perform weighting: at increasing dT, there is less data, so weigh
        % Binned_y based on the amount of data available:
        for iter = 1:size(Binned_y,1)
            Binned_y(iter,:) = Binned_y(iter,:).*JumpsPerdT(iter);
        end
    end
    
elseif ModelFit == 2
    %You want to fit to a CDF function
    %So first we must calculate the CDF from the finely binned PDF
    Binned_y_CDF = zeros(size(JumpProb,1), size(JumpProb,2));

    %Normalize the PDF
    for i=1:size(Binned_y_CDF,1)
        Binned_y_PDF(i,:) = y(i,:)./sum(y(i,:));
    end
    %calculate the CDF
    for i=1:size(Binned_y_CDF,1)
        for j=2:size(Binned_y_CDF,2)
            Binned_y_CDF(i,j) = sum(Binned_y_PDF(i,1:j));
        end
    end

    %Output the final variable
    Binned_y = Binned_y_CDF;
    
    if UseWeights == 1
        % Perform weighting: at increasing dT, there is less data, so weigh
        % Binned_y based on the amount of data available:
        for iter = 1:size(Binned_y,1)
            Binned_y(iter,:) = Binned_y(iter,:).*JumpsPerdT(iter);
        end    
    end
end

% ensure that the total sum of all fractions is one:
%   F_bound + F_free1 + F_free2 = 1
% an ad-hoc way of doing this is to check that F_bound + F_Free1 < 1
% if it is not, we add a very high cost to the output
if F_BOUND + F_FREE1 <= 1
    % OK, they are below one; so do not penalize:
    penalty_func = zeros(1, size(Binned_y,2));
elseif F_BOUND + F_FREE1 > 1
    penalty_func = 10^5*(1-F_BOUND - F_FREE1) .* ones(1, size(Binned_y,2));
end

% now append the custom cost function to the model output
Binned_y = vertcat(Binned_y, penalty_func);





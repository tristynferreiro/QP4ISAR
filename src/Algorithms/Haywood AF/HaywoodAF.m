function [AF_RA_HRRP] = HaywoodAF(RA_HRRP)
% Implements the Haywood autofocus algorithm
    % This is the phase compensation part of the Haywood Motion
    % Compensation algorithm. It uses the variance of power across all
    % range-aligned profiles to select a dominant scatterer (DS). The phase
    % of this DS is used to calculate the compensation factor to correct
    % all other scatterers.
    
    % A. Zyweck, PhD Thesis Appendix was used as a resource in
    % implementing this Haywood algorithm

    % Revision 1: Restrucutred the code for optimisation. Added additional
    % step comments

    %% Parameters
    num_range_bins = size(RA_HRRP,2);
    magnitude = abs(RA_HRRP);
    
    %% Step 1: Identify scatterer (range bin) that satifies the dominant scatterer criteria.
    
    % Criteria 1: power of scatterer>average power - Eq A.10 of Zyweck's appendix
    power_scatterer = sum(magnitude.^2,1);
    average_power_scatterer = mean(power_scatterer);
    
    % find possible scatterers using a threshold
    scalingFactor = 1;
    candidateScatterersIdx = find(power_scatterer>scalingFactor*average_power_scatterer); % array of matches
    
    % Calculate Eq A.8 of Zyweck's appendix - amplitude variance
    amplitude_variance = var(magnitude(:,candidateScatterersIdx),1); % 1xM matrix

    % Criteria 2: candidate with minimum variance is dominant scatterer (DS) - Eq A.9 of Zyweck's appendix
    [~,variance_DS_idx] = min(amplitude_variance); % returns index of match in candidateScatterersIdx
    DS_idx = candidateScatterersIdx(variance_DS_idx);
    
    %% Step 3: Calculate phase differences - Eq A.11 of Zyweck's appendix
    DS_phase_history = angle(RA_HRRP(:,DS_idx)); % N x 1 matrix
    
    %% Step 4: apply phase differences to each HRRP - Eq A.12 of Zyweck's appendix
    DS_complex_conjugate = exp(-1i*DS_phase_history);
    compensation_matrix = repmat(DS_complex_conjugate,1,num_range_bins);
    AF_RA_HRRP = RA_HRRP.* compensation_matrix; % Eq A.13 of Zyweck's appendix

end
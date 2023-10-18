function [AF_RA_HRRP] = HaywoodAF_v0(RA_HRRP)
% Implements the Haywood autofocus algorithm
    % This is the phase compensation part of the Haywood Motion
    % Compensation algorithm. It uses the variance of power across all
    % range-aligned profiles to select a dominant scatterer (DS). The phase
    % of this DS is used to calculate the compensation factor to correct
    % all other scatterers.
    %
    % A. Zyweck, PhD Thesis Appendix was used as a resource in
    % implementing this Haywood algorithm
    %
    % Version: v0

    %% Parameters
    num_range_bins = size(RA_HRRP,2);
    
    %% Step 1: Calculate Eq A.8 of Zyweck's appendix - amplitude variance
    amplitude_variance = var(abs(RA_HRRP), 1); % 1xM matrix
    
    %% Step 2: Identify scatterer (range bin) that satifies the dominant scatterer criteria.
    
    % Criteria 1: power of scatterer>average power - Eq A.10 of Zyweck's appendix
    power_scatterer = sum(abs(RA_HRRP).^2,1);
    average_power_scatterer = mean(power_scatterer);
    % find possible scatterers where power of scatterer > average power
    candidate_scatterers_idx = find(power_scatterer>average_power_scatterer); % array of matches
    
    % Criteria 2: candidate with minimum variance is dominant scatterer (DS) - Eq A.9 of Zyweck's appendix
    [~,variance_DS_idx] = min(amplitude_variance(candidate_scatterers_idx)); % returns index of match in candidateScatterersIdx
    DS_idx = candidate_scatterers_idx(variance_DS_idx);
    
    %% Step 3: Calculate phase differences - Eq A.11 of Zyweck's appendix
    DS_phase_history = angle(RA_HRRP(:,DS_idx)); % N x 1 matrix
    
    %% Step 4: apply phase differences to each HRRP - Eq A.12 of Zyweck's appendix
    DS_complex_conjugate = exp(-1i*DS_phase_history);
    compensation_matrix = repmat(DS_complex_conjugate,1,num_range_bins);
    AF_RA_HRRP = RA_HRRP.*compensation_matrix; % Eq A.13 of Zyweck's appendix
end
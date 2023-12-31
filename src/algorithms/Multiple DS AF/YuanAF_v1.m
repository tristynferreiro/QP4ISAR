function [AF_RA_HRRP] = YuanAF_v1(RA_HRRP)
% Implements the Yuan autofocus algorithm
    % It is a type of dominant scatterer algorithm (DSA) that uses multiple 
    % dominant scatterres (DS) to achieve phase correction. Candidate 
    % scatterers are selected based on a threshold value of 0.16 
    % (step 2). From these candidate scatterers a number of DSs are
    % selected (step 3). The phases of the DSs are used to calculate a 
    % compensation factor use to correct all other scatterers (step 4).    
    %
    % Chao Yuan and David P. Casasent, https://doi.org/10.1117/1.1425789
    % Note: A threshold scaling_factor was introduced to reduce effects of noise.
    %
    % Version: v1
    %
    % Revisions: 
    %   - Added if statement to set number of scatterers to highest number
    %     possible in cases where a smaller number of candidate scatterers
    %     are available.

    %% Variables
    num_range_bins = size(RA_HRRP,2);
    
    %% Step 1: Calculate mean and variance
    amplitude_mean = mean(abs(RA_HRRP),1);
    amplitude_variance = var(abs(RA_HRRP),1);
    
    %% Step 2:  find candidate scatterers: 
    % Identify candidate scatterers using Yuan's threshold
    dispersion = amplitude_variance./(amplitude_variance+amplitude_mean.^2);
    candidate_scatterers_idx  = find(dispersion< 0.16); % profile numbers
    
    % Plot Noisy Candidate Scatterers: HRRP profiles
    % noisy_HRRP = zeros(size(RA_HRRP,1), size(RA_HRRP,2));
    % noisy_HRRP(:,candidate_scatterers_idx) = RA_HRRP(:,candidate_scatterers_idx);
    % figure; imagesc(20*log10(abs(noisy_HRRP))); colormap('jet'); colorbar;
    % xlabel('Range (m)'); ylabel('Profile Number');
    % title('HRR Profiles: Selected scatterers (no noise filtering)');

    %% Step 3:  Choose smallest 11 (preferred) but can choose number between 6-18
    if(size(candidate_scatterers_idx,2)>11)
        num_scatterers = 11; % ideally 11 otherwise value in range 6-18
    else
        % set the number of scatterers as largest value value possible
        % given the number of candidate scatterers
        num_scatterers = size(candidate_scatterers_idx,2); 
    end

    [~,candidate_scatterers_idx_min] = mink(dispersion(candidate_scatterers_idx), num_scatterers);
    DS_idx = candidate_scatterers_idx(candidate_scatterers_idx_min); % get range bin numbers
    
    % Plot Noisy Dominant Scatterers: HRRP profiles
    % noisy_HRRP = zeros(size(RA_HRRP,1), size(RA_HRRP,2));
    % noisy_HRRP(:,DS_idx) = RA_HRRP(:,DS_idx);
    % figure; imagesc(20*log10(abs(noisy_HRRP))); colormap('jet'); colorbar;
    % xlabel('Range (m)'); ylabel('Profile Number');
    % title('HRR Profiles: Selected dominant scatterers (no noise filtering)');

    %% Step 4:  Determine constant phase shift for N pulses
    ref_bins = RA_HRRP(1,DS_idx); % reference profile
    product_vector = conj(ref_bins).* RA_HRRP(:,DS_idx);
    % calculate average phase difference across all range bins in each profile
    phase_shifts = angle(mean(product_vector,2)); 

    %% Step 5:  Apply phase shift
    compensation_angle = exp(-1i*phase_shifts);
    compensation_matrix = repmat(compensation_angle,1,num_range_bins);
    AF_RA_HRRP = RA_HRRP.*compensation_matrix;

end
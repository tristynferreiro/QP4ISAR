function [AF_RA_HRRP] = YuanAF(RA_HRRP)
% Implements the Yuan autofocus algorithm
    % This approach is based on the dominant scatterer algorithm (DSA). 

    num_range_bins = size(RA_HRRP,2);
    
    %% Step 1: Calculate mean and variance
    voltages = abs(RA_HRRP);
    amplitude_mean = mean(voltages,1);
    amplitude_variance = var(voltages,1);
    
    %% Step 2:  find candidate scatterers: 
    % Threshold to get target profiles
    no_noise_scatterers = find(amplitude_mean.^2>mean(amplitude_mean.^2));
    no_noise_scatterers = min(no_noise_scatterers):max(no_noise_scatterers);
    
    % Plot chosen HRRP profiles based on thresholding
    % test_HRRP = zeros(size(RA_HRRP,1), size(RA_HRRP,2));
    % test_HRRP(:,noNoiseScatterers) = RA_HRRP(:,noNoiseScatterers);
    % figure; imagesc(20*log10(abs(test_HRRP))); colormap('jet');
    % figure; imagesc(20*log10(abs(RA_HRRP))); colormap('jet');
    
    % Identify candidate scatterers using Yuan's threshold
    criteria = amplitude_variance./(amplitude_variance+amplitude_mean.^2);
    criteria_thresh = criteria(no_noise_scatterers)< 0.16;
    candidate_scatterers_idx = no_noise_scatterers(criteria_thresh);
    
    %% Step 3:  Choose smallest x scatterers
    if(size(candidate_scatterers_idx,2)>11)
        num_scatterers = 11; % ideally 11 otherwise value in range 6-18
    else
        % set the number of scatterers as largest value value possible
        % given the number of candidate scatterers
        num_scatterers = size(candidate_scatterers_idx,2); 
    end

    [~,candidate_scatterers_sorted] = sort(criteria(candidate_scatterers_idx));
    DS_idx = candidate_scatterers_idx(candidate_scatterers_sorted(1:num_scatterers));  % get range bin numbers
    
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
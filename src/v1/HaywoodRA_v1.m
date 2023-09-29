function [RA_HRRP] = HaywoodRA_v1(HRRP, ref_HRRP_num)
% Implements the Haywood range-alignment algorithm
    % This approach uses the correlation RA algorithm but incorporates
    % phase adjustment to achieve better range-alignment.


    %% Step 1 get the correlation values (the bin shifts): - Eq A.3 of Zyweck's appendix
    [~,shifts] = correlationRA(HRRP,ref_HRRP_num);

    %Plot stair case Function
    figure; plot(1:size(shifts,1),shifts)
    xlabel('Profile Number');
    ylabel('Number of bin shifts')
    title('Bin shifts per Range Profile');
    hold on

    %% Step 2: smooth the shift values so we polyfit to the lowest degree
    % polynomial
    [num_profiles,num_range_bins]=size(HRRP);
    profiles = 1:num_profiles;
    coefficients = polyfit(profiles,shifts,1);
    smoothed_coefficients = polyval(coefficients,profiles); % smoothed values

    % plot smoothed values
    plot(profiles,smoothed_coefficients,'-')
    hold off

    %% Step 3: Perform Haywood's Range alignment
    m = 0:1:num_range_bins-1;
    % calculate phase correction angle - Eq A.6 of Zyweck's appendix
    phase = exp(-1i*(2*pi*smoothed_coefficients'.*m)/num_range_bins);
    % Apply correction angle to profiles - - Eq A.5 of Zyweck's appendix
    RA_HRRP = ifft(phase.*fft(HRRP,num_range_bins,2),num_range_bins,2);
end
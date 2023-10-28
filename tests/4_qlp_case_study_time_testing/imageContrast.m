function [image_contrast] = imageContrast(image)
% Calculate the image contrast (IC) of the input image. 
    % Uses Eq3 of M. Martorella and F. Berizzi, "Time windowing for highly 
    % focused isar image reconstruction," IEEE Transactions on Aerospace 
    % and Electronic Systems, vol. 41, no. 3, pp. 992–1007, 2005

num_elements = numel(image);
%% Step 1: Calculate intensity = amplitude^2
intensity = abs(image).^2;

%% Step 2: Calculate the a mean average power
pixel_mean = sum(intensity,"all")/num_elements;

%% Step 3: Calculate the IC
% difference and square
numerator = sum((intensity - pixel_mean).^2,"all")/num_elements;
% calculate contrast value
image_contrast = sqrt(numerator)/pixel_mean;
end
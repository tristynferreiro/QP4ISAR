function [image_contrast] = imageContrast(image)
% Calculate the image contrast (IC) as described in 
% M. Martorella and F. Berizzi, “Time windowing for highly focused isar image reconstruction,” IEEE
% Transactions on Aerospace and Electronic Systems, vol. 41, no. 3, pp. 992–1007, 2005

numElements = numel(image);
%% Step 1: Calculate intensity = amplitude^2
intensity = abs(image).^2;

%% Step 2: Calculate the a mean average power
pixel_mean = sum(intensity,"all")/numElements;

%% Step 3: Calculate the IC
% difference and square
numerator = sum((intensity - pixel_mean).^2,"all")/numElements;
% calculate contrast
image_contrast = sqrt(numerator)/pixel_mean;
end
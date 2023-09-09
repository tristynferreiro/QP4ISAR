function [contrast] = imageContrast(imageArray)
% Calculate the image contrast

numElements = numel(imageArray);

% calculate intensity=amplitude^2
intensity = real(imageArray).^2;
% calculate the mean
pixel_mean = sum(intensity,"all")/numElements;
% difference and square
numerator = sum((intensity - pixel_mean).^2,"all")/numElements;
% calculate contrast
contrast = sqrt(numerator)/(sum(intensity,"all")/numElements);
end
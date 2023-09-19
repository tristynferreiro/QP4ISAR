function [contrast] = imageContrast(image)
% Calculate the image contrast

numElements = numel(image);

% calculate intensity=amplitude^2
intensity = abs(image).^2;
% calculate the mean
pixel_mean = sum(intensity,"all")/numElements;
% difference and square
numerator = sum((intensity - pixel_mean).^2,"all")/numElements;
% calculate contrast
contrast = sqrt(numerator)/(sum(intensity,"all")/numElements);
end
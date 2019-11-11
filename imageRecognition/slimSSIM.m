%% SSIM function designed to work with greyscale 300 x 300 images
% 

function val = slimSSIM(im1,im2)

[x,y] = meshgrid(1:3,1:3);
window = exp(-(...
    ((x-ceil(length(x)/2)).^2)./(2*1.5^2) +...
    ((y-ceil(length(y)/2)).^2)./(2*1.5^2)...
    ));
window = window/sum(window,'all');
window = rot90(window,2);

% c(2) = 58.5225;
% c(1) = 6.5025;
c(2) = (0.05*255)^2;
c(1) = (0.02*255)^2;

avg1 = conv2(im1,window,'valid');
avg2 = conv2(im2,window,'valid');
avg11 = avg1.^2;
avg22 = avg2.^2;
avg12 = avg1.*avg2;
s11 = conv2(im1.^2,window,'valid') - avg11;
s22 = conv2(im2.^2,window,'valid') - avg22;
s12 = conv2(im1.*im2,window,'valid') - avg12;

mat = ((2*avg12 + c(1)).*(2*s12 + c(2)))./...
    ((avg11 + avg22 + c(1)).*(s11 + s22 + c(2)));

val = mean(mat(:));

end
function [signal] = unbias_image(img)

    [m,n] = size(img) ;
    % simply fit a line at boundary
    bias = img(:,1) ;
    [P] = polyfit((1:1:m)',bias,1) ;
    slope = P(1) ;
    intercept = P(2) ;
    bias = slope*(1:1:m)'+intercept ;
    background = bias*ones(1,n) ;
    signal = img - background ;

end


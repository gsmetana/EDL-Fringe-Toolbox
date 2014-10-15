function f = fexpand(g,mm,nn)
    % expand f to [m n]
    [m,n] = size(g) ;
    f0 = g ;
    % generate a larger matrix with size [mm nn]
    f = zeros(mm,nn) ;
    % copy original data
    f(1:m,1:n) = f0 ;
end


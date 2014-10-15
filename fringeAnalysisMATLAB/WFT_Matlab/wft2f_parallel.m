clear all
close all
clc

tic

pathname = [pwd '/Images/'] ;
files = dir([pathname '/' '*.tif']) ;
l = length(files) ;

sigmax = 10 ;
sigmay = 10 ;
wxl = -.5 ;
wxh = 0 ;
wxi = 0.05;
wyl = -.5 ;
wyh = .5 ;
wyi = 0.05 ;

% PARAMETERS FOR LOCATING SPHERE IN EACH IMAGE
% pix_m = 26667 ;  % PIX/METER
% dt = 1e-4 ;      % TIME STEP
% Vsphere = 2.4 ;  % SPHERE SPEED
% c0 = [463 124] ; % SPHERE CENTER AT IMAGE 34
% r = 87 ;         % SPHERE RADIUS IN PIXELS

for k = 1:l
    num = k ;
    t = num*dt ;
    g = double(imread([pathname files(num).name])) ;
    bg = imopen(g,strel('disk',7)) ;
    img = g-bg ;
    img = imadjust(img) ;
    a = mean2(img) ;
    f = img-a ;
    % SPHERE CENTER OF EACH IMAGE
    % c = [c0(1) round(c0(2)+Vsphere*(num-34)*dt*pix_m)] ; 
    % f = MaskSphere(f,r,c) ;
    sx = round(3*sigmax) ;
    sy = round(3*sigmay) ;  
    
    % image size
    [m,n] = size(f) ;
    mm = m+2*sx ;
    nn = n+2*sy ;
    f = fexpand(f,mm,nn) ;
    Ff = fft2(f) ;
    [y,x] = meshgrid(-sy:1:sy,-sx:1:sx+1) ;
    w0 = exp(-x.*x/sigmax/sigmax/2-y.*y/sigmay/sigmay/2) ;
    w0 = w0/sum(sum(w0.*w0)).^(1/2) ;
       
    yN = wyl:wyi:wyh ;
    xN = wxl:wxi:wxh ;
    leny = length(yN) ;
    lenx = length(xN) ;
    lent = leny*lenx*1.00 ;
    A = zeros(1,m*n) ;

    for i = 1:lenx
        Z = zeros(1,m*n) ;
        perc = [num2str(i/lenx*100) '% complete'] ;
        disp(perc)
        parfor j=1:leny
            wxt = xN(i) ;
            wyt = yN(j) ;
            w = w0.*exp(1i*wxt*x+1i*wyt*y) ;
            w = fexpand(w,mm,nn) ;
            Fw = fft2(w) ;
            sf = ifft2(Ff.*Fw) ;
            sf = sf(1+sx:m+sx,1+sy:n+sy) ;
            [M,N] = size(sf) ;
            sf_reshape = reshape(sf,[1,M*N]) ;
            Z = [Z;sf_reshape] ;          
        end
        [M,N] = size(Z) ;
        Sf = Z(2:M,1:N) ;
        [maxSf] = max(Sf,[],1) ;
        A = max([A;maxSf],[],1) ;
    end
    
    display('End of parallel computing')
    display([num2str(num) ' of ' num2str(l) ' images complete']) 

    f = f(1:m,1:n) ;
    maxSf = reshape(A,[m,n]) ;
    g.r = abs(maxSf) ;
    g.phase = angle(maxSf) ;
    g.filtered = g.r.*exp(g.phase*sqrt(-1)) ;
    % g.filtered = MaskSphere(g.filtered,r,c) ;
    unwrapped = unwrapping_qg_trim(g.filtered) ;
    phase = unbias_image(unwrapped) ;
    % phase = MaskSphere(phase,r,c) ;
    
    xtick = 0:200:n ;
    ytick = 0:200:m ;
    figure;
    subplot(2,1,1) 
    imagesc(f);
    colormap gray ;
    title('Disturbed Image $g(x,y)$','FontSize',16) ;
    set(gca,'FontSize',12,'XTick',xtick,'YTick',ytick)
    colorbar
    axis image ;
    subplot(2,1,2) 
    imagesc(phase);
    title('Signal','FontSize',16) ;
    set(gca,'FontSize',12,'XTick',xtick,'YTick',ytick)
    axis image ;
    colorbar
    colormap jet ;
    colormap(flipud(colormap))
    caxis([min(min(phase)) 0])
    set(gcf,'PaperSize',[4 7],'PaperPosition',[0 0.1 4 7]) ;    
    saveas(gcf,['Results/Images/Unwrapped/UnwrappedPhase_' num2str(num) '.pdf'],'pdf')

    dlmwrite(['Results/Data/Amplitude/' num2str(num) '.txt'],abs(g.filtered)) ;
    dlmwrite(['Results/Data/Phase/' num2str(num) '.txt'],angle(g.filtered)) ;
    dlmwrite(['Results/Data/Unwrapped/' num2str(num) '.txt'],phase) ;
end

toc
                    
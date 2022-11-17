clc
clear all
close all
 

PperM = 1389.52/(7*12e-3);

fid = fopen('C:\Users\banv20\OneDrive - University of Bath\Bath - Postgraduate\Post FYP\X_Drive\FEA_results.txt', 'wt');
for n = 1:114
    %% ORIGINAL IMAGE
    
    left = 'C:\Users\banv20\OneDrive - University of Bath\Bath - Postgraduate\Post FYP\X_Drive\FEA Images\';
    right = '.png';
    join_name = join([left, num2str(n), right]);
    carrot = imread(join_name);
    
     
    %% ISOLATING COLOUR CHANNELS AND GETTING THERE GREYSCALE EQUIVALENCE, AND APPLY OTSU THRESHOLDING:
     
    % BLUE CHANNEL
    blue = carrot(:,:,3);
    levelB = graythresh(blue);
    otsublue = levelB*255;
    binarycarrotB = im2bw(blue,levelB);
    
     
    % RED CHANNEL
    red = carrot(:,:,1);
    levelR = graythresh(red);
    otsured = levelR*255;
    binarycarrotR = im2bw(red,levelR);
     
    %% IMAGE PROCESSING
     
    % Inverting binary
    binarycarrotBi = 1 - binarycarrotB;
    binarycarrot = binarycarrotBi.*binarycarrotR;
    imshow(binarycarrot);
    SE = strel('diamond',7);
    CI = imerode(binarycarrot,SE);
    CI = imdilate(CI,SE);
     
    %% KEEPING THE LARGEST BLOB 

    x = imbinarize(CI);
    [a bc]=bwlabel(x);
    b=[];
    temp=0;
    gk=1;
    [f g]=size(a);
    for i=1:bc
        for wk=1:f
            for tk=1:g
                if(a(wk,tk)==gk)
                    temp=temp+1;
                end
            end
        end
        b=[b temp];
        temp=0;
        gk=gk+1;
    end
    [m n]=max(b);
     
    output=[];
    for i=1:f
        for j=1:g
            if(a(i,j)==n)
                output(i,j)=1;
            else
                output(i,j)=0;
            end
        end
    end
    subplot(1,2,1);
    imshow(x);
    title('Original Image');
    subplot(1,2,2);
    imshow(output);
    title('Keeping largest object in the binary image');
     
     
    %% GET THE UPPER EDGE OF BLOB
     
    [sizey, sizex] = size(output);
    xpix = 1:sizex;
    ypix = zeros(1,sizex);
    for n = 1:sizex
        if sum(output(:,n)) > 1
            ypix(n) = find(output(:,n) == 1, 1, 'first');
        end
    end
    %correct plot for the image size in meters
    xval = (xpix/PperM);
    yval = ((sizey - ypix)/PperM);
    %Check the plot
    figure(11)
    plot(xval,yval,'r')
    xlim([0 (sizex/(PperM))])
    ylim([0 (sizey/(PperM))])
    %find the region that the carrot occupies
     
    datmin = find(ypix ~= 0, 1, 'first');
    datmax = find(ypix ~= 0, 1, 'last');
    ytrue = yval(datmin:datmax);
    xtrue = xval(datmin:datmax)-xval(datmin);
     
    %Check the plot
    figure(11)
    plot(xtrue(1:6),ytrue(1:6),'g','linewidth',2)
    hold on
    plot(xtrue(6:end),ytrue(6:end),'r','linewidth',2)
    xlim([0 (sizex/(PperM))])
    ylim([0 (sizey/(PperM))])
    hold off
     
    %% Fit Circle
     
    x=xtrue(:); y=ytrue(:);
    a=[x y ones(size(x))]\[-(x.^2+y.^2)];
    xc = -.5*a(1);
    yc = -.5*a(2);
    R_auto  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
    p = nsidedpoly(1000, 'Center', [xc yc], 'Radius', R_auto);
    hold on
    plot(p, 'FaceColor', 'b')
    hold off
    
    fprintf(fid, '%s', num2str(R_auto));
    fprintf(fid, '\n');

end
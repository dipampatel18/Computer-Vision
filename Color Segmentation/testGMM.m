% Initialization
clear;
close all;
clc;

k = 7; % Number of the gaussians
[mu,covar,phi] = trainGMM_trial1(k);
f = measureDepth();
threshold = 0.94;
pdf = [];
% Training
cd test_images;
imgfiles = dir('*.jpg');
num_imgs = length(imgfiles);

for i = 1:num_imgs %23 images
    img = imgfiles(i).name;
    curr_img = imread(img);
    %curr_img = imgaussfilt(curr_img);
    curr_img = imsharpen(curr_img);
    
    curr_img = im2double(curr_img); % Converts the intensity image I to double precision, rescaling the data if necessary
    [l,b,d] = size(curr_img);
    
    % Seperate the RGB channels 
    img_channel_R = curr_img(:,:,1);
    img_channel_G = curr_img(:,:,2);
    img_channel_B = curr_img(:,:,3);

    % Matrix to vector conversion
    img_channel_R = img_channel_R(:);
    img_channel_G = img_channel_G(:);
    img_channel_B = img_channel_B(:);
    
    % Genearting the input data, 307200(640X480)X3
    X = [img_channel_R,img_channel_G,img_channel_B];
    
    % Calculating the likelihood using formula
%    [len,width] = size(X);
%    y = zeros(len,1);
%     for i = 1 : len
%         temp2 = (X(i,:)-mu);
%         temp2 = temp2*covar_inv;
%         temp2 = temp2*((X(i,:)-mu)');
%         temp2 = exp((-0.5)*temp2);
%         y(i)=temp1*temp2;
%     end

    % Computing the likelihood and the posterior for the input images
    for i = 1 : k
        curr_covar = cell2mat(covar(1,i));
        pdf = mvnpdf(X,mu(i,:),curr_covar);
        pdf = phi(i)*pdf;
        y = reshape(pdf,l,b);
        final_pdf = cat(3,y);
    end 
      prob_pixel = sum(final_pdf,3);     

     % Genearting a black and white image to locate the ball
    BW = im2bw(curr_img);
    indexW = (prob_pixel >= threshold);
    indexB = (prob_pixel < threshold);
    BW(indexW) = 1;
    BW(indexB) = 0;
%     BW = bwmorph(BW,'close');
%     BW = bwmorph(BW,'fill');

    % Removing pixels based on the number of connected pixels
    LB = 50; % Lower bound
    UB = 1890; % Upper bound
    BW = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));
    
    % Estimate Center and Radii of Circular Objects and Plot Circles
    stats = regionprops('table',BW,'Centroid','MajorAxisLength','MinorAxisLength');
    centers = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2;
    area = pi*(radii.^2); % Computing the area of the ball
    distance_temp = max(f(area)); % Distance estimation based on the model
    if isempty(distance_temp)
        distance = 0;
    else
        distance = distance_temp;
    end
    position = [350 320];
    curr_img = insertText(curr_img,position,distance,'FontSize',18,'TextColor','black');
    
    figure, imshow(BW);
    figure, imshow(curr_img);
    hold on
    viscircles(centers,radii);
    hold off
%    imshowpair(curr_img, BW, 'montage')    
end
cd ..
% Plotting the GMM ellipsoids 
figure;
for i = 1:k
        plotGMM((mu(i,:))',cell2mat(covar(1,i)));
        hold on;
end



% Initialization
clear;
close all;
clc;

%mu = [];
%covar = [];
%area = zeros(23,2);
threshold = 0.94; % The connfidence score
prob_orange = 0.2; % Prior, orange or not orange
[mu,covar] = trainSingleGaussian();
f = measureDepth(); 

% Values to be used in likelihood estimation
covar_inv = inv(covar);
covar_det = det(covar);
temp_var1 = (2*pi)^3;
temp_var1 = sqrt(temp_var1*covar_det);
temp_var1 = 1/temp_var1;

cd test_images;
imgfiles = dir('*.jpg');
num_imgs = length(imgfiles);

for i = 1:num_imgs %23 images
    
    img = imgfiles(i).name;
    curr_img = imread(img);
    curr_img = imgaussfilt(curr_img);
    curr_img = im2double(curr_img); % Converts the intensity image I to double precision, rescaling the data if necessary
    
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
    
    % Computing the likelihood using formula
    [len,~] = size(X);
    y = zeros(len,1);
    for i = 1 : len
        temp_var2 = (X(i,:)-mu);
        temp_var2 = temp_var2*covar_inv;
        temp_var2 = temp_var2*((X(i,:)-mu)');
        temp_var2 = exp((-0.5)*temp_var2);
        y(i)=temp_var1*temp_var2;
    end
    
     [l,b,~] = size(curr_img);
     prob1 = reshape(y,l,b);
     
     % Computing the likelihood using inbuilt matlab function
     prob2 = mvnpdf(X,mu,covar);
     prob2 = reshape(prob2,l,b);
     
    % Computing the Posterior
    prob_pixel = prob2 * prob_orange;
    
    % Genearting a black and white image to locate the ball
    BW = im2bw(curr_img);
    indexW = (prob_pixel >= threshold);
    indexB = (prob_pixel < threshold);
    BW(indexW) = 1;
    BW(indexB) = 0;
%     BW = bwmorph(BW,'close');
%     BW = bwmorph(BW,'fill');

    % Removing pixels based on the number of connected pixels
    LB = 117; % Lower bound
    UB = 1890; % Upper bound
    BW = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));
    
    % Estimate Center and Radii of Circular Objects and Plot Circles
    stats = regionprops('table',BW,'Centroid','MajorAxisLength','MinorAxisLength');
    centers = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2;
    area = pi*(radii.^2); % area of the ball 
    distance = max(f(area)); % Estimating the distance of the ball given the regression model
    
    % Distance of the ball marked in the image
    position = [300 320];
    curr_img = insertText(curr_img,position,distance,'FontSize',10,'TextColor','black');
    
    figure;
    imshow(curr_img);
    hold on
    viscircles(centers,radii);
    hold off
%     imshowpair(curr_img, BW_improved, 'montage');
     
end
cd ..


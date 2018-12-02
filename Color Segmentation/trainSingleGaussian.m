function [mu,covar] = trainSingleGaussian()
cd Cropped_Imgs;
imgfiles = dir('*.jpg');
num_imgs = length(imgfiles);

mu = [];
covar = [];

for i = 1:num_imgs %23 images
    img = imgfiles(i).name;
    curr_img = imread(img);
    curr_img = imgaussfilt(curr_img);
    curr_img = im2double(curr_img); % Converts the intensity image I to double precision, rescaling the data if necessary
    
    % Seperate the RGB channels 
    cropped_img_R = curr_img(:,:,1);
    cropped_img_G = curr_img(:,:,2);
    cropped_img_B = curr_img(:,:,3);
    
    % Matrix to vector conversion
    cropped_img_R = cropped_img_R(:);
    cropped_img_G = cropped_img_G(:);
    cropped_img_B = cropped_img_B(:);
    
    % Computing mu for each channel
    mu_cropped_img_R = mean(cropped_img_R);
    mu_cropped_img_G = mean(cropped_img_G);
    mu_cropped_img_B = mean(cropped_img_B);
    mu_temp = [mu_cropped_img_R;mu_cropped_img_G;mu_cropped_img_B]; % 3X1, for a single image
    mu = [mu mu_temp]; % 3X23, for 23 images
    
    % Computing the covariance matrix
    % covariance = [std_R^2  cov_RG  cov_RB;
    %               cov_RG std_G^2   cov_BG;
    %               cov_RB cov_BG  std_B^2];
    % Computing standard deviation for each channel
    std_R = std(cropped_img_R);
    std_G = std(cropped_img_G);
    std_B = std(cropped_img_B);
    % Computing covariances between channels
    cov_RG = cov(cropped_img_R,cropped_img_G);
    %cov_GR = cov(cropped_img_G,cropped_img_R);
    cov_RB = cov(cropped_img_R,cropped_img_B);
    %cov_BR = cov(cropped_img_B,cropped_img_R);
    cov_BG = cov(cropped_img_B,cropped_img_G);
    %cov_GB = cov(cropped_img_G,cropped_img_B);
    cov_temp = [std_R^2  cov_RG(1,2)  cov_RB(1,2);cov_RG(2,1) std_G^2   cov_BG(1,2);cov_RB(2,1) cov_BG(2,1)  std_B^2]; %3X3
    covar = cat(3,covar,cov_temp); %3X3X23
end
cd ..
% Computing the average mean and covariance for 23 images
mu = mean(mu,2);
mu = mu';
covar = mean(covar,3);
end
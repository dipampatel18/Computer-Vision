function [mu_GMM,covar_GMM,phi_GMM] = trainGMM_trial1(num_Gaussian)
threshold = 0.001;
X = [];
covar = [];
% length_img = 0;
% width_img = 0;
X = data();

% Parameter initialization to random values
m = size(X, 1); % Number of data points.
k = num_Gaussian;  % Number of Gaussians
n = 3;  % Dimension of the data, 3 in our case

% Randomly selected k data points to serve as the initial means.
indx = randperm(m);
mu = X(indx(1:k), :);

% Overall covariance of the dataset as the initial variance for each gaussian.
for j = 1 : k
    covar{j} = cov(X);
end

% weightage of the gaussians(prior probabilities for gaussians) 
phi = ones(1, k) * (1 / k);

% Expectation Maximization
% Weightage of each datapoint to the gaussians
alpha = [];
fprintf('Training in progress');
for iterations = 1:200 % Change the value for number of iterations 
    %%%%%%%%Expectation step%%%%%%%%%%%%%%%%%%%%
    % Likelihood values for each data point for every gaussian.
    prob = []; 
     % For each Gaussian
    for j = 1 : k
        % Computing the likelihood for all data points for cluster 'j'.
        prob(:, j) = mvnpdf(X, mu(j, :), covar{j});
    end
    
    % Multiply each prob value by the gaussian weight(prior) to get the posterior
    %  prob  [m  x  k]  phi  [1  x  k]  prob_pixel  [m  x  k]    
    % Computing the gaussian weights
    alpha = (prob.*phi) ./ sum((prob.*phi),2);
    
    %%%%%%%%%%Maximization step%%%%%%%%%%%%%%%%%
    prevMu = mu; % previous mean to evaluate the convergence criterion
    % Average weight/prior.
    temp = mean(alpha, 1);
    for (j = 1 : k) % compute mu,covar and wt/prior for each cluster j
        % Restricting the value of phi between 0.1 and 0.21 for better
        % results
        if (temp(j)<0.1 || temp(j)>0.21)
            phi(j) = 0.2;
        else
            phi(j) = temp(j);
        end
        
        % Weighted average to get the new mean
        mu(j, :) = ((alpha(:, j))' * X)/ sum(alpha(:, j));
        
        covar_temp = zeros(n, n);
        Xm = X - mu(j,:);
        % Compute covariance for each pixel and take summation
        for (i = 1 : m)
            covar_temp = covar_temp + (alpha(i, j) .* (Xm(i, :)' * Xm(i, :)));
        end
        % Weighted average to get covariance matrix
        covar{j} = covar_temp / sum(alpha(:, j));
    end
     % Convergence criterion
    diff = mu-prevMu;
    if (diff>=threshold)
        fprintf('Loop ended after iterations: %d', iterations);
        break
    end
    
end
fprintf('Loop ended after iterations: %d', 200);
mu_GMM = mu;
covar_GMM = covar;
phi_GMM = phi;

% % Multivariant normal probability density function
% function mult_pdf = my_mvnpdf(multi_X,multi_mu,multi_covar)
%     % Values to be used in likelihood estimation
%     covar_inv = inv(multi_covar);
%     covar_det = det(multi_covar);
%     temp_var = (2*pi)^3;
%     temp_var = sqrt(temp_var*covar_det);
%     temp_var = 1/temp_var;
%     % Computing the likelihood using formula
%     [len,~] = size(multi_X);
%     y = zeros(len,1);
%     for i = 1 : len
%         temp_var2 = (multi_X(i,:)-multi_mu);
%         temp_var2 = temp_var2*covar_inv;
%         temp_var2 = temp_var2*(temp_var2');
%         temp_var2 = exp((-0.5)*temp_var2);
%         y(i)=temp_var*temp_var2;
%     end
%     mult_pdf = reshape(y,length_img,width_img);
% end

% Generating training data
function X = data()
x = [];
cd Cropped_Imgs;
imgfiles = dir('*.jpg');
num_imgs = length(imgfiles);
for i = 1:num_imgs %23 training images
    img = imgfiles(i).name;
    curr_img = imread(img);
    curr_img = imgaussfilt(curr_img);
    curr_img = im2double(curr_img); % Converts the intensity image I to double precision, rescaling the data if necessary
%     [length_img,width_img,~] = size(curr_img);
    % Seperate the RGB channels 
    curr_img_R = curr_img(:,:,1);
    curr_img_G = curr_img(:,:,2);
    curr_img_B = curr_img(:,:,3);
    
%    image = cat(3, curr_img_R, curr_img_G, curr_img_B);
%    imshowpair(image, maskImage, 'montage')
    
    % Taking transpose
    curr_img_R = curr_img_R';
    curr_img_G = curr_img_G';
    curr_img_B = curr_img_B';
    
    % Matrix to vector conversion
    cropped_img_R = curr_img_R(:);
    cropped_img_G = curr_img_G(:);
    cropped_img_B = curr_img_B(:);
    
    x = [x;cropped_img_R,cropped_img_G,cropped_img_B];
end
X = x;
end
cd ..
end

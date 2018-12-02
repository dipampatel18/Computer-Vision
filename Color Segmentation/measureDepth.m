function f = measureDepth()
x = []; % Input features, area in pixels
y = []; % Distance of the ball in centimeter

cd Cropped_Imgs;
imgfiles = dir('*.jpg');
num_imgs = length(imgfiles);

for i = 1:num_imgs %23 images
    % Generating the training data set
    img = imgfiles(i).name;
    curr_img = imread(img);
    [l,b,~] = size(curr_img);
    area = l*b;
    x = [x;area];
    val = str2double(extractBefore(img,4));
    y = [y;val];
end
% figure; % open a new figure window
% plot(x, y, 'rx', 'MarkerSize', 10); % Plot the data
% ylabel('Distance of the ball in cm'); 
% xlabel('Region area in pixels'); 
f=fit(x,y,'poly3','Normalize','on','Robust','Bisquare'); % Fits a quadratic function to the data set
%plot(f,x,y);
cd ..
end


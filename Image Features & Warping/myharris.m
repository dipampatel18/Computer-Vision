%% CMSC 426: COMPUTER VISION

% Name: DIPAM PATEL
% UID: 115809833
% Homework-2: Image Features and Warping

function corner_response = myharris(I, window_size, corner_thresh)

clear all;
close all;
clc;


%% Image-1

%% Reading, Denoising and Converting the Image to Gray

Img = imread ('1.jpg');

Img = medfilt3(Img);

I = rgb2gray(Img);

% Initialization of Variables

k = 0.04;                       % The value of k ranges from 0.04 to 0.06

corner_thresh = 10000000;        % Threshold value selected after trial and error from the R values

sigma = 1;                      % Sigma of Gaussian

window_size = 1;                % Window under consideration having a standard size of 1

order = (2*window_size + 1);    % To find the localMaxima

%% Applying the Prewitt Filter

[Gx, Gy] = meshgrid(-window_size:window_size, -window_size:window_size);

% The Prewitt Operator has been used here where Gx = [1 0 -1; 1 0 -1; 1 0 -1]
% and Gy = Gx'

%% Computing the x and y Derivatives of the Image

Ix = conv2(Gx, I);      % Convoluting to control the size of the output and for direct comparison with the input
Iy = conv2(Gy, I);      
                        
%% Computing Products of Derivatives at Every Pixel
                        
Ix2 = Ix .* Ix;         % Multiplying the Derivates to every element
Iy2 = Iy .* Iy;         
Ixy = Ix .* Iy;         

%% Compute the sums of the products of derivatives at each pixel

G = [1 1 1; 1 1 1; 1 1 1];      % Taking I(3x3) matrix for convoluting over the product of derivates at every pixel.

Sx2 = conv2(G, Ix2);
Sy2 = conv2(G, Iy2);
Sxy = conv2(G, Ixy);

%% Calculations for Every Pixel

[r, c] = size(I);

for x = 1:r
   
    for y = 1:c
       
% Forming a Matrix at Every Pixel
        
        H = [Sx2(x,y) Sxy(x,y);
             Sxy(x,y) Sy2(x,y)];
       
% Computing the Response of the Detector at each Pixel
        
        R = det(H) - k*(trace(H) ^ 2);
       
% Threshold on value of R

        if (R > corner_thresh)
            im(x,y) = R;        % Forming a matrix having only values higher than the threshold
        
        else
            im(x,y) = 0;
        end

   end
   
end

%% Computing the Local Maxima of the given Pixels

localMaxima = ordfilt2(im, order^2, ones(order));   % Out of the 9 Pixels present inside the window of consideration 
                                                    % ordfilt2(x, 9, ones(3*3))

harrisPoints = (im >= localMaxima) & (im > corner_thresh);  % Finding the Harris Points based on the value of Local Maxima and Threshold value
[rows, cols] = find(harrisPoints);          % Finding a vector containing the linear indices of each nonzero element in the array

%% Plotting the Harris Corners 

figure, imshow(Img);
hold on;

plot(cols, rows, 'r+');
hold off;

%% Plotting the Heatmap

figure,
imshow(im);



%% Image-2

clear all;
clc;

%% Reading, Denoising and Converting the Image to Gray

Img = imread ('2.jpg');

Img = medfilt3(Img);

I = rgb2gray(Img);

% Initialization of Variables

k = 0.04;                       % The value of k ranges from 0.04 to 0.06

corner_thresh = 10000000;        % Threshold value selected after trial and error from the R values

sigma = 1;                      % Sigma of Gaussian

window_size = 1;                % Window under consideration having a standard size of 1

order = (2*window_size + 1);    % To find the localMaxima

%% Applying the Prewitt Filter

[Gx, Gy] = meshgrid(-window_size:window_size, -window_size:window_size);

% The Prewitt Operator has been used here where Gx = [1 0 -1; 1 0 -1; 1 0 -1]
% and Gy = Gx'

%% Computing the x and y Derivatives of the Image

Ix = conv2(Gx, I);      % Convoluting to control the size of the output and for direct comparison with the input
Iy = conv2(Gy, I);      
                        
%% Computing Products of Derivatives at Every Pixel
                        
Ix2 = Ix .* Ix;         % Multiplying the Derivates to every element
Iy2 = Iy .* Iy;         
Ixy = Ix .* Iy;         

%% Compute the sums of the products of derivatives at each pixel

G = [1 1 1; 1 1 1; 1 1 1];      % Taking I(3x3) matrix for convoluting over the product of derivates at every pixel.

Sx2 = conv2(G, Ix2);
Sy2 = conv2(G, Iy2);
Sxy = conv2(G, Ixy);

%% Calculations for Every Pixel

[r, c] = size(I);

for x = 1:r
   
    for y = 1:c
       
% Forming a Matrix at Every Pixel
        
        H = [Sx2(x,y) Sxy(x,y);
             Sxy(x,y) Sy2(x,y)];
       
% Computing the Response of the Detector at each Pixel
        
        R = det(H) - k*(trace(H) ^ 2);
       
% Threshold on value of R

        if (R > corner_thresh)
            im(x,y) = R;        % Forming a matrix having only values higher than the threshold
        
        else
            im(x,y) = 0;
        end

   end
   
end

%% Computing the Local Maxima of the given Pixels

localMaxima = ordfilt2(im, order^2, ones(order));   % Out of the 9 Pixels present inside the window of consideration 
                                                    % ordfilt2(x, 9, ones(3*3))

harrisPoints = (im >= localMaxima) & (im > corner_thresh);  % Finding the Harris Points based on the value of Local Maxima and Threshold value
[rows, cols] = find(harrisPoints);          % Finding a vector containing the linear indices of each nonzero element in the array

%% Plotting the Harris Corners 

figure, imshow(Img);
hold on;

plot(cols, rows, 'g*');
hold off;

%% Plotting the Heatmap

figure,
imshow(im);


end

%% References

% 1. http://www.cse.psu.edu/~rtc12/CSE486/lecture06.pdf
% 2. http://alumni.media.mit.edu/~maov/classes/vision09/lect/09_Image_Filtering_Edge_Detection_09.pdf
% 3. https://cmsc426.github.io/pano-prereq/
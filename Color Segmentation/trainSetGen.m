clear;
close all;
clc;

cd train_images;
imgfiles = dir('*.jpg');
for i = 1:length(imgfiles)
    % Read the given training images
    img = imgfiles(i).name;
    %imshow(imread(img));
    curr_img = imread(img);
    
    % Displays the image in a figure window and creates an interactive Crop Image tool 
    % associated with the image
    cropped_img = imcrop(curr_img);
    
    % Save the cropped training images showing just the orange ball in Cropped_Imgs folder
    baseFileName = img; % Same file name as the read image
    fullFileName = fullfile('D:\Robotics\CMSC426\P1\train_images\Cropped_Imgs',baseFileName);
    imwrite(cropped_img,fullFileName);    
    
end
cd ..
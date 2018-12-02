%% CMSC 426: COMPUTER VISION

% Name: DIPAM PATEL
% UID: 115809833
% Homework-1: Linear Least Squares

%%
%% Plotting Eigenvalues and Eigenvectors
%%

clear all
close all
clc

N = 200;                        % Sample Size

%% Data-1

% Loading the .mat Files into MATLAB 

load('data1.mat');

x1 = pts(1,:);                  % Assigning first row of the matrix to x 
y1 = pts(2,:);                  % and second row to y from the .mat file

figure;
plot(x1, y1, '*');                              % Plotting to visualize the data
title('Data-1: Eigenvalues and Eigenvectors');  % Title for the scatter plot
hold on;

mean_x1 = mean(x1);             % Finding the Mean of the Data
mean_y1 = mean(y1);

temp_x1 = 0;                    % Initialization
temp_y1 = 0;
temp_xy1 = 0;

%% Calculating the Covariance Matrix, Eigenvector and Eigenvalue 

% Using the Formula of Calculating the Covariance as mentioned in Theory

for i = 1:N
    
    x_temp1 = x1(i) - mean_x1;  
    y_temp1 = y1(i) - mean_y1;
    
    tempx1 = x_temp1 * x_temp1;
    tempy1 = y_temp1 * y_temp1;
    tempxy1 = x_temp1 * y_temp1;
    
    temp_x1 = temp_x1 + tempx1;
    temp_y1 = temp_y1 + tempy1;
    temp_xy1 = temp_xy1 + tempxy1;
    
end

% Formation of Covariance Matrix

EV1 = [temp_x1 temp_xy1;        % Covariance Matrix = [xx xy; xy yy]
      temp_xy1 temp_y1];

% Eigenvector and Eigenvalue

[V1, D1] = eig(EV1);            % Diagonal matrix D of eigenvalues and
                                % matrix V whose columns are the corresponding
                                % right eigenvectors, so that A*V = V*D.


%% In the Direction of Maximum Eigenvalue

Maximum = max(max(D1));         % Inner 'max' gives the maximum column of
                                % the matrix and outer 'max' gives the
                                % maximum from that column

[xmax1, ymax1] = find(D1 == Maximum);   % 'find' returns the row and column
                                        % subscripts of each nonzero element
                                        % in the given array using any of the input
                                        % arguments in previous syntaxes.

UV1 = V1(:, ymax1);             % Taking the maximum column value from the eigenvector
slope1 = UV1(2,1)/UV1(1,1);     % Using the value of UV, calculating the
                                % direction of the vector to be plotted

xrangemax1 = -100:100;              % Scaling the Magnitude of Vector for
                                    % better visualization of the axes
yrangemax1 = slope1 * xrangemax1;   

plot(xrangemax1, yrangemax1)
hold on;

%% In the Direction of Minimum Eigenvalue

Minimum = min(min(D1));         % Inner 'max' gives the maximum column of
                                % the matrix and outer 'max' gives the
                                % maximum from that column

[xmin1, ymin1] = find(D1 == Minimum);   % 'find' returns the row and column
                                        % subscripts of each nonzero element
                                        % in the given array using any of the input
                                        % arguments in previous syntaxes.

UV1 = V1(:,ymin1);              % Taking the maximum column value from the eigenvector
slope1 = UV1(2,1)/UV1(1,1);     % Using the value of UV, calculating the
                                % direction of the vector to be plotted

xrangemin1 = -5:5;                  % Scaling the Magnitude of Vector for
                                    % better visualization of the axes
yrangemin1 = slope1 * xrangemin1;

plot(xrangemin1, yrangemin1)
hold off;


%%
%% Data-2
%%

% Loading the .mat Files into MATLAB 

load('data2.mat');

x2 = pts(1,:);               
y2 = pts(2,:);               

figure;
plot(x2, y2, '*');                              
title('Data-2: Eigenvalues and Eigenvectors');  
hold on;

mean_x2 = mean(x2);          
mean_y2 = mean(y2);

temp_x2 = 0;                 
temp_y2 = 0;
temp_xy2 = 0;

%% Calculating the Covariance Matrix, Eigenvector and Eigenvalue 

for i = 1:N
    
    x_temp2 = x2(i) - mean_x2;
    y_temp2 = y2(i) - mean_y2;
    
    tempx2 = x_temp2 * x_temp2;
    tempy2 = y_temp2 * y_temp2;
    tempxy2 = x_temp2 * y_temp2;
    
    temp_x2 = temp_x2 + tempx2;
    temp_y2 = temp_y2 + tempy2;
    temp_xy2 = temp_xy2 + tempxy2;
    
end

% Formation of Covariance Matrix

EV2 = [temp_x2 temp_xy2;      
      temp_xy2 temp_y2];

% Eigenvector and Eigenvalue

[V2, D2] = eig(EV2);           

%% In the Direction of Maximum Eigenvalue

Maximum = max(max(D2));

[xmax2, ymax2] = find(D2 == Maximum);

UV2 = V2(:, ymax2);
slope2 = UV2(2,1)/UV2(1,1);

xrangemax2 = -100:100;
yrangemax2 = slope2 * xrangemax2;

plot(xrangemax2, yrangemax2)
hold on;

%% In the Direction of Minimum Eigenvalue

Minimum = min(min(D2));

[xmin2, ymin2] = find(D2 == Minimum);

UV2 = V2(:,ymin2);
slope2 = UV2(2,1)/UV2(1,1);

xrangemin2 = -5:5;
yrangemin2 = slope2 * xrangemin2;

plot(xrangemin2, yrangemin2)
hold off;


%%
%% Data-3
%%

% Loading the .mat Files into MATLAB 

load('data3.mat');

x3 = pts(1,:);       
y3 = pts(2,:);       

figure;
plot(x3, y3, '*');   
title('Data-3: Eigenvalues and Eigenvectors');
hold on;

mean_x3 = mean(x3);  
mean_y3 = mean(y3);

temp_x3 = 0;         
temp_y3 = 0;
temp_xy3 = 0;

%% Calculating the Covariance Matrix, Eigenvector and Eigenvalue 

for i = 1:N
    
    x_temp3 = x3(i) - mean_x3;
    y_temp3 = y3(i) - mean_y3;
    
    tempx3 = x_temp3 * x_temp3;
    tempy3 = y_temp3 * y_temp3;
    tempxy3 = x_temp3 * y_temp3;
    
    temp_x3 = temp_x3 + tempx3;
    temp_y3 = temp_y3 + tempy3;
    temp_xy3 = temp_xy3 + tempxy3;
    
end

% Formation of Covariance Matrix

EV3 = [temp_x3 temp_xy3;
      temp_xy3 temp_y3];

% Eigenvector and Eigenvalue

[V3, D3] = eig(EV3);    

%% In the Direction of Maximum Eigenvalue

Maximum = max(max(D3));

[xmax3, ymax3] = find(D3 == Maximum);

UV3 = V3(:, ymax3);
slope3 = UV3(2,1)/UV3(1,1);

xrangemax3 = -100:100;
yrangemax3 = slope3 * xrangemax3;

plot(xrangemax3, yrangemax3);
hold on;

%% In the Direction of Minimum Eigenvalue

Minimum = min(min(D3));

[xmin3, ymin3] = find(D3 == Minimum);

UV3 = V3(:,ymin3);
slope3 = UV3(2,1)/UV2(1,1);

xrangemin3 = -5:5;
yrangemin3 = slope3 * xrangemin3;

plot(xrangemin3, yrangemin3)
hold off;


%% References

% 1. https://cmsc426.github.io/math-tutorial/
% 2. http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
% 3. https://en.wikipedia.org/wiki/Random_sample_consensus
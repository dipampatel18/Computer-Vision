%% CMSC 426: COMPUTER VISION

% Name: DIPAM PATEL
% UID: 115809833
% Homework-1: Linear Least Squares

%%
%% Line Fitting & Outliers Rejection
%%

clear all
close all
clc

%% Data-1

% Loading the '.mat' File into MATLAB

load('data1.mat');

x1 = pts(1,:);              % Assigning first row of the matrix to x 
y1 = pts(2,:);              % and second row to y


%% Linear Least Square

figure;
subplot(3,3,2);
plot(x1, y1, '*');
title('Data-1: Linear Least Square');
hold on;

trans_x1 = x1';             % Transpose of x and y and converting them into
trans_y1 = y1';             % column vectors

trans_x1 = [ones(size(x1', 1), 1), trans_x1];
                            % For ease of calculation and to match the 
                            % matrix dimensions, the column matrix is
                            % converted into 200*2 with first column of
                            % ones only
                                                
theta1 = (inv(trans_x1' * trans_x1)) * (trans_x1' * trans_y1);
                            % Using the formula for Calculating the 'theta'
                            % parameter as mentioned in the Theory

plot(trans_x1(:, 2), trans_x1*theta1, 'o');
hold on;


%% Outlier Rejection using Regularization

% Case-1: Lambda = 10:

subplot(3,3,4);
plot(x1, y1, '*');
title('Data-1: Outlier Rejection using Regularization, L = 10');
hold on;

lambda1 = 10;               % Varying the value of Lambda to see its effect
                            % on the plot
                            
theta_reg1 = (inv((trans_x1' * trans_x1) + (lambda1 * eye(size(trans_x1, 2)))) * (trans_x1' * trans_y1));
                            % Using the formula for Calculating the 'theta'
                            % parameter as mentioned in the Theory

plot(trans_x1(:, 2), trans_x1 * theta_reg1, '+');


% Case-2: Lambda = 100:

subplot(3,3,5);
plot(x1, y1, '*');
title('Data-1: Outlier Rejection using Regularization, L = 100');
hold on;

lambda1 = 100;              % Varying the value of Lambda to see its effect
                            % on the plot
                            
theta_reg1 = (inv((trans_x1' * trans_x1) + (lambda1 * eye(size(trans_x1, 2)))) * (trans_x1' * trans_y1));
                            % Using the formula for Calculating the 'theta'
                            % parameter as mentioned in the Theory

plot(trans_x1(:, 2), trans_x1 * theta_reg1, '+');


% Case-3: Lambda = 1000:

subplot(3,3,6);
plot(x1, y1, '*');
title('Data-1: Outlier Rejection using Regularization, L = 1000');
hold on;

lambda1 = 1000;             % Varying the value of Lambda to see its effect
                            % on the plot
                            
theta_reg1 = (inv((trans_x1' * trans_x1) + (lambda1 * eye(size(trans_x1, 2)))) * (trans_x1' * trans_y1));
                            % Using the formula for Calculating the 'theta'
                            % parameter as mentioned in the Theory

plot(trans_x1(:, 2), trans_x1 * theta_reg1, '+');

%% Outlier Rejection using RANSAC

subplot(3,3,8);
plot(x1, y1, '*');
title('Data-1: Outlier Rejection using RANSAC');
hold on;

comp_count1 = 0;            % Initializing Counter

for i=1:1000
    
    rand_int_11 = randi(size(trans_x1, 1)); % Randomly selecting the points
                                            % for iterations
    rand_int_12 = randi(size(trans_x1, 1));

    a = (y1(:, rand_int_12) - y1(:, rand_int_11));
    b = (x1(:, rand_int_11) - x1(:, rand_int_12));
    c = (x1(:, rand_int_11) * (y1(:, rand_int_11) - y1(:, rand_int_12)))+(y1(:, rand_int_11)*(x1(:,rand_int_12)-x1(:,rand_int_11)));
                                            % Setting up the variables of
                                            % the equation of line
    dist1 = abs(((a*x1) + (b*y1) + c)/sqrt((a^2) + (b^2)));
                                            % Calculating the perpendicular
                                            % distance from the point to line
    
    thres1 = 5;             % Setting up the threshold by trial and error
    
    dist_thres1 = dist1 < thres1;
    count_temp1 = sum(sum(dist_thres1));

    if count_temp1 > comp_count1
        
        comp_count1 = count_temp1;
        
        index_11 = rand_int_11;
        index_12 = rand_int_12;
        
        b_final1 = b;
        a_final1 = a;
        c_final1 = c;
        
    end

end

slope_ransac1 = (y1(:, index_12) - y1(:, index_11))/(x1(:, index_12) - x1(:, index_11));
                            % Calculating the slope of the RANSAC line
                            % using the index variable as column for x and y
                            
x_ransac_1 = -100:100;      
y_ransac_1 = [-c_final1 - (a_final1 * x_ransac_1)]/(b_final1);
                            % Determining the value of y by keeping a range
                            % for x
                            
plot(x_ransac_1, y_ransac_1,'*');   % Plotting the RANSAC line
hold off;


%%
%% Data-2
%%

% Loading the '.mat' File into MATLAB

load('data2.mat');

x2 = pts(1,:);             
y2 = pts(2,:);             


%% Linear Least Square

figure;
subplot(3,3,2);
plot(x2, y2, '*');
title('Data-2: Linear Least Square');
hold on;

trans_x2 = x2';            
trans_y2 = y2';            

trans_x2 = [ones(size(x2', 1), 1), trans_x2];
theta2 = (inv(trans_x2' * trans_x2)) * (trans_x2' * trans_y2);
                           
plot(trans_x2(:, 2), trans_x2 * theta2, 'o');
hold on;


%% Outlier Rejection using Regularization

% Case-1: Lambda = 10:

subplot(3,3,4);
plot(x2, y2, '*');
title('Data-2: Outlier Rejection using Regularization, L = 10');
hold on;

lambda2 = 10;
theta_reg2 = (inv((trans_x2' * trans_x2) + (lambda2 * eye(size(trans_x2,2))))*(trans_x2'*trans_y2));
plot(trans_x2(:, 2), trans_x2 * theta_reg2, '+');


% Case-2: Lambda = 100:

subplot(3,3,5);
plot(x2, y2, '*');
title('Data-2: Outlier Rejection using Regularization, L = 100');
hold on;

lambda2 = 100;
theta_reg2 = (inv((trans_x2' * trans_x2) + (lambda2 * eye(size(trans_x2,2))))*(trans_x2'*trans_y2));
plot(trans_x2(:, 2), trans_x2 * theta_reg2, '+');


% Case-3: Lambda = 1000:

subplot(3,3,6);
plot(x2, y2, '*');
title('Data-2: Outlier Rejection using Regularization, L = 1000');
hold on;

lambda2 = 1000;
theta_reg2 = (inv((trans_x2' * trans_x2) + (lambda2 * eye(size(trans_x2,2))))*(trans_x2'*trans_y2));
plot(trans_x2(:, 2), trans_x2 * theta_reg2, '+');


%% Outlier Rejection using RANSAC

subplot(3,3,8);
plot(x2, y2, '*');
title('Data-2: Outlier Rejection using RANSAC');
hold on;

comp_count2 = 0;

for i=1:1000
    
    rand_int_21 = randi(size(trans_x2, 1));
    rand_int_22 = randi(size(trans_x2, 1));


    a = (y2(:, rand_int_22) - y2(:, rand_int_21));
    b = (x2(:, rand_int_21) - x2(:, rand_int_22));
    c = (x2(:, rand_int_21) * (y2(:, rand_int_21) - y2(:, rand_int_22))) + (y2(:, rand_int_21) * (x2(:,rand_int_22) - x2(:,rand_int_21)));
    
    dist2 = abs(((a*x2) + (b*y2) + c)/sqrt((a^2) + (b^2)));
    thres2 = 5;
    
    dist_thres2 = dist1 < thres2;
    count_temp2 = sum(sum(dist_thres2));

    if count_temp2 > comp_count2
        
        comp_count2 = count_temp2;
        
        index_21 = rand_int_21;
        index_22 = rand_int_22;
        
        b_final2 = b;
        a_final2 = a;
        c_final2 = c;
        
    end

end

slope_ransac2 = (y2(:, index_22) - y2(:, index_21))/(x2(:, index_22) - x2(:, index_21));

x_ransac_2 = -100:100;
y_ransac_2 = [-c_final2 - (a_final2 * x_ransac_2)]/(b_final2);

plot(x_ransac_2, y_ransac_2,'*');
hold off;


%%
%% Data-3
%%

% Loading the '.mat' File into MATLAB

load('data3.mat');

x3 = pts(1,:);               
y3 = pts(2,:);               


%% Linear Least Square

figure;
subplot(3,3,2);
plot(x3, y3, '*');
title('Data-3: Linear Least Square');
hold on;

trans_x3 = x3';              
trans_y3 = y3';              

trans_x3 = [ones(size(x3', 1), 1), trans_x3];
theta3 = (inv(trans_x3' * trans_x3)) * (trans_x3' * trans_y3);
                           
plot(trans_x3(:, 2), trans_x3 * theta3, 'o');
hold on;


%% Outlier Rejection using Regularization

% Case-1: Lambda = 10:

subplot(3,3,4);
plot(x3, y3, '*');
title('Data-3: Outlier Rejection using Regularization, L = 10');
hold on;

lambda3 = 10;
theta_reg3 = (inv((trans_x3' * trans_x3) + (lambda3 * eye(size(trans_x3, 2))))*(trans_x3' * trans_y3));
plot(trans_x3(:, 2), trans_x3 * theta_reg3, '+');

% Case-1: Lambda = 100:

subplot(3,3,5);
plot(x3, y3, '*');
title('Data-3: Outlier Rejection using Regularization, L = 100');
hold on;

lambda3 = 100;
theta_reg3 = (inv((trans_x3' * trans_x3) + (lambda3 * eye(size(trans_x3, 2))))*(trans_x3' * trans_y3));
plot(trans_x3(:, 2), trans_x3 * theta_reg3, '+');

% Case-1: Lambda = 1000:

subplot(3,3,6);
plot(x3, y3, '*');
title('Data-3: Outlier Rejection using Regularization, L = 1000');
hold on;

lambda3 = 1000;
theta_reg3 = (inv((trans_x3' * trans_x3) + (lambda3 * eye(size(trans_x3, 2))))*(trans_x3' * trans_y3));
plot(trans_x3(:, 2), trans_x3 * theta_reg3, '+');

%% Outlier Rejection using RANSAC

subplot(3,3,8);
plot(x3, y3, '*');
title('Data-1: Outlier Rejection using RANSAC');
hold on;

comp_count3 = 0;

for i=1:1000
    
    rand_int_31 = randi(size(trans_x3, 1));
    rand_int_32 = randi(size(trans_x3, 1));


    a = (y3(:, rand_int_32) - y3(:, rand_int_31));
    b = (x3(:, rand_int_31) - x3(:, rand_int_32));
    c = (x3(:, rand_int_31) * (y3(:, rand_int_31) - y3(:, rand_int_32)))+(y3(:, rand_int_31)*(x3(:,rand_int_32)-x3(:,rand_int_31)));
    
    dist3 = abs(((a*x3) + (b*y3) + c)/sqrt((a^2) + (b^2)));
    thres3 = 5;
    
    dist_thres3 = dist3 < thres3;
    count_temp3 = sum(sum(dist_thres3));

    if count_temp3 > comp_count3
        
        comp_count3 = count_temp3;
        
        index_31 = rand_int_31;
        index_32 = rand_int_32;
        
        b_final3 = b;
        a_final3 = a;
        c_final3 = c;
        
    end

end

slope_ransac3 = (y3(:, index_32) - y3(:, index_31))/(x3(:, index_32) - x3(:, index_31));

x_ransac_3 = -100:100;
y_ransac_3 = [-c_final3 - (a_final3 * x_ransac_3)]/(b_final3);

plot(x_ransac_3, y_ransac_3,'*');
hold off;

%% References

% 1. https://cmsc426.github.io/math-tutorial/
% 2. http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
% 3. https://en.wikipedia.org/wiki/Random_sample_consensus
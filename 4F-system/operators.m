%% ????????????
clear
load imdemos circuit
%I = rgb2gray(imread('ex2','jpg')); 
I=imread('a','bmp');
I=imresize(I,[300 300]);
I = double(I);
imshow(I, []);
%h = [0 -1 0;-1 4 -1; 0 -1 0];%laplacian
%h = [-1 -2 -1;0 0 0; 1 2 1];%sobel
h = [0 1 0;-1 0 1; 0 -1 0];%one derivation
%h = [-1 -1 -1;-1 8 -1; -1 -1 -1];
J = conv2(I, h, 'same');
K = I - J;
figure, imshow(J, []);


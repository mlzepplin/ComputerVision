% Demo of "Region Based Active Contours"
%
% Example:
% seg_demo
%
% Coded by: Shawn Lankton (www.shawnlankton.com)
% Updated by: Rishabh Gupta (mlzepplin)

imageFile = '\\Client\C$\Users\rishabh\Desktop\texture.jpg';
I = imread(imageFile);
m = zeros(size(I,1),size(I,2));          %-- create initial mask
m(75:175,65:165) = 1;

I = imresize(I,.5);  %-- make image smaller 
m = imresize(m,.5);  %     for fast computation

figure(3);

subplot(2,2,1); imshow(I); title('Input Image');
subplot(2,2,2); imshow(m); title('Initialization');
subplot(2,2,3); title('Segmentation');

seg = region_seg(I, m, 2000); %-- Run segmentation

subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');



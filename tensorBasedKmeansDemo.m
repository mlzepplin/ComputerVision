clear;
%img = imread('\\Client\H$\Desktop\texture.jpg');
[filename,pathname]=uigetfile('*.jpg','Select the image');
img = imread(fullfile(pathname,filename), 'jpg');
form = 1;
clusters = tensorBasedKmeans(img,3,form);
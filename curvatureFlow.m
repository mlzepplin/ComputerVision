clear;

%curve given as a png image of 1 pixel width
[filename,pathname]=uigetfile('*.png','Select the image');
img = imread(fullfile(pathname,filename), 'png');

%another way to input file
%imageFile = '\\Client\C$\Users\rishabh\Desktop\rish.png';
%img = imread(imageFile);

imgGray = rgb2gray(img);
iw = imgGray<255;
B = bwboundaries(iw);
B=B(1);
B = B{1,1}';

xvalues = B(1,:);
yvalues = B(2,:);


%CIRCLE/ ELLIPSE
% theta = 0:pi/1000:2*pi;
% theta = theta(1,1:size(theta,2)-1);
% xvalues = 100*cos(theta);
% yvalues = 100*sin(theta);


figure(1);
plot(xvalues,yvalues,'b');
hold on

   xx = xvalues;
   yy = yvalues;

    for i=1:10000 
        xp = deriv(xvalues);
        yp = deriv(yvalues);
        xpp = deriv(xp);
        ypp = deriv(yp);
        k = ((xpp.*yp)-(ypp.*xp))./(sqrt(1+ xp.^2 + yp.^2)).^3;
        nx = yp./(sqrt(1+ xp.^2 + yp.^2));
        ny = -xp./(sqrt(1+ xp.^2 + yp.^2));
        cxt = nx.*(k);
        cyt = ny.*(k);
        xvalues = xvalues + 0.5*cxt;
        yvalues = yvalues + 0.5*cyt;
        if mod(i,1000)== 0
            plot(xvalues,yvalues,'g');
            hold on
            drawnow;
        end
       
        i
    end
plot(xvalues,yvalues,'r');
    
%taking derivative by finite differences method
function x=deriv(a)
x=zeros(size(a));
length = size(a,2);
        for j=1:length
                if j==1
                    x(1,j) = (a(1,j+1)-a(1,length))/2;
                elseif j==length
                    x(1,j) = (a(1,1)-a(1,j-1))/2;
                else
                    x(1,j) = (a(1,j+1)-a(1,j-1))/2; 
                end
        end

end





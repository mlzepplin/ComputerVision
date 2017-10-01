clear
imageFile = '\\Client\C$\Users\rishabh\Desktop\marston.jpeg';
img = imread(imageFile);
imgGray = double(rgb2gray(img));
figure('Name','Input image');imagesc(imgGray);
%size(XY)
%size(XYPrime)
XY = Get2DPoints(imageFile, 4);
XYPrime = Get2DPoints(imageFile, 4); 

A = computeA(XY,XYPrime,4);
[U,S,V] = svd(A);
[m,n] = size(V);
h = V(:,n);
H = reshape(h,[3,3]);
H=H';
Hinv = inv(H);
[row,col] = size(imgGray);
INew = zeros(row,col);

XYPrime
for i=1:col
    for j=1:row
        imgOldPoint =Hinv*[i,j,1]'; %i,j are new image coordinates
        
        %inverse mapping and using bilinear interpolation
        ximgOldPoint = imgOldPoint(1,1)/imgOldPoint(3,1);
        yimgOldPoint = imgOldPoint(2,1)/imgOldPoint(3,1);
        xOld = floor(ximgOldPoint);
        yOld = floor(yimgOldPoint);
        a = abs(ximgOldPoint - xOld);
        b = abs(yimgOldPoint - yOld);
        
        if xOld+1<col && xOld>0 && yOld+1<row && yOld>0  
        intensityOld = (1-b)*(1-a)*imgGray(yOld,xOld) +  b*(1-a)*imgGray(yOld+1,xOld) +a*b*imgGray(yOld+1,xOld+1)+(1-b)*a*imgGray(yOld,xOld+1);
        else
            intensityOld = 0;
        end
        INew(j,i) = intensityOld;
    end
end
figure('Name','Output Image');imagesc(INew);

function A = computeA(XY,XYPrime,numPoints)
    if size(XY)~=size(XYPrime)
        return  
    else
        X = XY(1,:);
        Y = XY(2,:);

        XPrime = XYPrime(1,:);
        YPrime = XYPrime(2,:);
        A = computeAi(X(1,1),Y(1,1),XPrime(1,1), YPrime(1,1));
        for i=2:numPoints
           ATemp = computeAi(X(1,i),Y(1,i),XPrime(1,i), YPrime(1,i));
           A = vertcat(A,ATemp);
        end
    end
end

function Ai= computeAi(x,y,xPrime,yPrime)
    Ai = [-x,-y,-1,0,0,0,x*xPrime,y*xPrime,xPrime;
           0,0,0,-x,-y,-1,x*yPrime,y*yPrime,yPrime];
    return;   
end

function XY2D=Get2DPoints(ImageFileName, NumberOfPoints)
XY2D=[];
[Img, Col]=imread(ImageFileName,'jpg'); %%% assuming the imagefile is a jpg file.
image(Img); drawnow; hold on;
for i=1:NumberOfPoints
[x, y]=ginput(1);
v=[x;y];
plot(x, y, 'r*');
XY2D=[XY2D v];
end
return;
end

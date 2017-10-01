clear
imageFile = '\\Client\C$\Users\rishabh\Desktop\q3';


n=6;
XYZReal = [[0;1;3],[0;1;2],[0;1;1],[0;1;0],[1;0;0],[2;0;0]];

XYImage= Get2DPoints(imageFile, n);
%let C be the caliberation matrix,computing A from Image = CReal, such that
%A[column vector C] =0, so that I can use SVD of A and get the vector c
%from it's null space. C = [c1 c2 c3 c4;c5 c6 c7 c8; c9 c10 c11 c12], and
%c=[c1,c2,c3,...,c11,c12]'

%A is the matrix that converts MX = b problem into AX=0 homogeneous
%equation
A = computeA(XYZReal,XYImage,n);
[U,S,V] = svd(A);
[m,n] = size(V);
[min,index] = min(diag(S(1:12,1:12)));

m = V(1:12, index);
M = reshape(m,[4,3]);

%SO M IS THE TOTAL CALIBERATION MATRIX
M=M'
M = M/(M(3,4));

%taking out other parameters
%M = (KR|KRT) = (E|ET), K=internal parameter matrix, R = rotaition matrix,
%T = Translation matrix

E = M(1:3,1:3);

%using RQ decomposition of E to get K and R, E=AB, A - upper triangular
%matrix, corresponds to K and B corresponds to R

%APPLYING RQ DECOMPOSITION ON E IS SAME AS APPLYING QR DECOMPOSITION ON E'
[A,B] = qr(E);

K = A;
R = B;

%and finally trnslation vector
T = inv(E)*M(:,4);

%Hence we've gotten all the matrices, and hence all the coefficients of caliberation matrix M
K
R
T

function A = computeA(XYZReal,XYImage,numPoints)

        X = XYZReal(1,:);
        Y = XYZReal(2,:);
        Z = XYZReal(3,:);
        
        XImage = XYImage(1,:);
        YImage = XYImage(2,:);
        A = computeAi(X(1,1),Y(1,1),Z(1,1),XImage(1,1), YImage(1,1));
        for i=2:numPoints
           ATemp = computeAi(X(1,i),Y(1,i),Z(1,i),XImage(1,i), YImage(1,i));
           A = vertcat(A,ATemp);
        end
    
end


function Ai= computeAi(xReal,yReal, zReal, xImage,yImage)
    Ai = [xReal,yReal,zReal,1,0,0,0,0,-xImage*xReal,-xImage*yReal,-xImage*zReal,-xImage;
           0,0,0,0,xReal,yReal,zReal,1,-yImage*xReal,-yImage*yReal,-yImage*zReal,-yImage];
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

% Tensor based Region Based Active Contour Segmentation
%
% seg = region_seg(I,init_mask,max_its,alpha,display)
%
% Inputs: I           2D image
%         init_mask   Initialization (1 = foreground, 0 = bg)
%         max_its     Number of iterations to run segmentation for
%         alpha       (optional) Weight of smoothing term
%                       higer = smoother.  default = 0.2
%         display     (optional) displays intermediate outputs
%                       default = true
%
% Outputs: seg        Final segmentation mask (1=fg, 0=bg)
%
% Description: This code implements the paper: "Active Contours Without
% Edges" By Chan Vese. This is a nice way to segment images whose
% foregrounds and backgrounds are statistically different and homogeneous.
%
%Update: on top of chan-vesse, implementing the paper: "An Affine Invariant Tensor Dissimilarity 
%Measure and its Applications to Tensor-valued Image Segmentation"
%
% Example:
% img = imread('tire.tif');
% m = zeros(size(img));
% m(33:33+117,44:44+128) = 1;
% seg = region_seg(img,m,500);
%
% Coded by: Shawn Lankton (www.shawnlankton.com)
% Updated by Rishabh Gupta (mlzepplin)
%------------------------------------------------------------------------
function seg = region_seg(I,init_mask,max_its,alpha,display)
  
   %NOTE----- IMPORTANT----
   %Set this form variable to 1 for case one of question
   %and set it to 2 for the J-divergance case
   form = 1;

  %-- default value for parameter alpha is .1
  if(~exist('alpha','var')) 
    alpha = .2; 
  end
  %-- default behavior is to display intermediate outputs
  if(~exist('display','var'))
    display = true;
  end
  %-- ensures image is 2D double matrix
  I = im2graydouble(I);  
%   g = fspecial('Gaussian');
%   I = imfilter(I, g);
%   
  
   sigma = 0.5;
    hsize = 3;
    h = fspecial('gaussian',hsize,sigma);
    I = conv2(I,h,'same');
  % computing the structure tensors
%   % getting Ix and Iy
    h = fspecial('sobel')./8;
    Iy = imfilter(I,h); 

    h = h';
    Ix = imfilter(I,h);
    
    %creating tensor field of the image
    STMatrix = zeros(size(I,1),size(I,2),2,2);
    for i=1:size(I,1)
        for j=1:size(I,2)
            STMatrix(i,j,1,1) = Ix(i,j)^2;
            STMatrix(i,j,1,2) = Ix(i,j)*Iy(i,j);
            STMatrix(i,j,2,1) = Ix(i,j)*Iy(i,j);
            STMatrix(i,j,2,2) = Iy(i,j)^2;
        end
    end
    
    %appying averaging by gaussian to take care of the 
    %sigularities
    %making 4 tensor images
    tensorI1 = zeros(size(I,1),size(I,2));
    tensorI2 = zeros(size(I,1),size(I,2));
    tensorI3 = zeros(size(I,1),size(I,2));
    tensorI4 = zeros(size(I,1),size(I,2));
    for i=1:size(I,1)
        for j=1:size(I,2)
            tensorI1(i,j) = STMatrix(i,j,1,1);
            tensorI2(i,j) = STMatrix(i,j,1,2);
            tensorI3(i,j) = STMatrix(i,j,2,1);
            tensorI4(i,j) =STMatrix(i,j,2,2);
        end
    end
    sigma = 0.5;
    hsize = 4;
    h = fspecial('gaussian',hsize,sigma);
    tensorI1 = conv2(tensorI1,h,'same');
    tensorI2 = conv2(tensorI2,h,'same');
    tensorI3 = conv2(tensorI3,h,'same');
    tensorI4 = conv2(tensorI4,h,'same');
    
    %now updating the structure tensor matrix with this
     for i=1:size(I,1)
        for j=1:size(I,2)
            STMatrix(i,j,1,1) = tensorI1(i,j);
            STMatrix(i,j,1,2) = tensorI2(i,j);
            STMatrix(i,j,2,1) = tensorI3(i,j);
            STMatrix(i,j,2,2) = tensorI4(i,j);
        end
     end
    
    %Now STMatrix is the neighbourhood smoothed structure 
    %tensor field of I
%   
  %-- Create a signed distance map (SDF) from mask
  phi = mask2phi(init_mask);
 
  
  %--main loop
  for its = 1:max_its   % Note: no automatic convergence test

    [idx,idy] = find(phi <=1.2 & phi >= -1.2);  %get the curve's narrow band
    temp = find(phi <= 1.2 & phi >= -1.2);
    
    %-- find interior and exterior mean
    [uptsx,uptsy] = find(phi<=0);                 % interior points
    [vptsx,vptsy] = find(phi>0);                  % exterior points
    
    u = meanFunc(STMatrix,uptsx,uptsy,form); %interior mean
    v = meanFunc(STMatrix,vptsx,vptsy,form); %exterior mean
            
    %UPDATING THIS FOR STRUCTURE TENSORS INSTEAD OF INTENSITY
    %sum of distance of inside points
    F = distanceSum(STMatrix,idx,idy,u,v,form);

     
    curvature = get_curvature(phi,temp);
   
    dphidt = F./max(abs(F)) + alpha*curvature;  % gradient descent to minimize energy
    
    %-- maintain the CFL condition
    dt = .45/(max(dphidt)+eps);
        
    %-- evolve the curve
    phi(temp) = phi(temp) + dt.*dphidt;

    %-- Keep SDF smooth
    phi = sussman(phi, .5);

    %-- intermediate output
    if((display>0)&&(mod(its,20) == 0)) 
      showCurveAndPhi(I,phi,its);  
    end
  end
  
  %-- final output
  if(display)
    showCurveAndPhi(I,phi,its);
  end
  
  %-- make mask from SDF
  seg = phi<=0; %-- Get mask from levelset

  
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------
function u = meanFunc(st,ptsx,ptsy,form)

  if form == 1  
      u = zeros(2,2);
      for i=1:2
            for j=1:2
                 for a=1:length(ptsx)
                 u(i,j) = u(i,j) + st(ptsx(a,1),ptsy(a,1),i,j)/(length(ptsx)+eps); % interior mean
                 end
            end
       end     
  elseif form == 2
      %computing A and B
      A = zeros(2,2);
      B = zeros(2,2);
      for i=1:length(ptsx)
          A(:,:)= A(:,:) + reshape(st(ptsx(i),ptsy(i),:,:),2,2);
          B(:,:)=B(:,:) + inv(reshape(st(ptsx(i),ptsy(i),:,:),2,2));
      end
      u = sqrtm(inv(B))*(sqrtm(sqrtm(B)*A*sqrtm(B)))*sqrtm(inv(B));
      
  end
  
  
function F = distanceSum(st,idx,idy,u,v,form)
     F = zeros(length(idx),1);
   if form == 1   
    
     for b=1:length(idx)
        for i=1:2
            for j=1:2
            F(b,1) = F(b,1) + (st(idx(b,1),idy(b,1),i,j)-u(i,j))^2-(st(idx(b,1),idy(b,1),i,j)-v(i,j))^2; % exterior mean
            
            end
        end
     end
   elseif form == 2
         for i=1:length(idx)
          F(i,1)= getJDist(reshape(st(idx(i),idy(i),:,:),2,2),u)-getJDist(reshape(st(idx(i),idy(i),:,:),2,2),v);
  
         end 
   end

    function d2 = getJDist(tensor,u)
      d2 = 1/2*trace(inv(tensor)*u +inv(u)*tensor)-2;
      
      
%-- Displays the image with curve superimposed
function showCurveAndPhi(I, phi, i)
  imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
  contour(phi, [0 0], 'g','LineWidth',4);
  contour(phi, [0 0], 'k','LineWidth',2);
  hold off; title([num2str(i) ' Iterations']); drawnow;
  
%-- converts a mask to a SDF
function phi = mask2phi(init_a)
  phi=bwdist(init_a)-bwdist(1-init_a);% im2double(init_a);
  
%-- compute curvature along SDF
function curvature = get_curvature(phi,idx)
    [dimy, dimx] = size(phi);        
    [y x] = ind2sub([dimy,dimx],idx);  % get subscripts

    %-- get subscripts of neighbors
    ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

    %-- bounds checking  
    ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;    

    %-- get indexes for 8 neighbors
    idup = sub2ind(size(phi),yp1,x);    
    iddn = sub2ind(size(phi),ym1,x);
    idlt = sub2ind(size(phi),y,xm1);
    idrt = sub2ind(size(phi),y,xp1);
    idul = sub2ind(size(phi),yp1,xm1);
    idur = sub2ind(size(phi),yp1,xp1);
    iddl = sub2ind(size(phi),ym1,xm1);
    iddr = sub2ind(size(phi),ym1,xp1);
    
    %-- get central derivatives of SDF at x,y
    phi_x  = -phi(idlt)+phi(idrt);
    phi_y  = -phi(iddn)+phi(idup);
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
    phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
             +0.25*phi(iddr)+0.25*phi(idul);
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;
    
    %-- compute curvature (Kappa)
    curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
              (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
  
%-- Converts image to one channel (grayscale) double
function img = im2graydouble(img)    
  [dimy, dimx, c] = size(img);
  if(isfloat(img)) % image is a double
    if(c==3) 
      img = rgb2gray(uint8(img)); 
    end
  else           % image is a int
    if(c==3) 
      img = rgb2gray(img); 
    end
    img = double(img);
  end

%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
  % forward/backward differences
  a = D - shiftR(D); % backward
  b = shiftL(D) - D; % forward
  c = D - shiftD(D); % backward
  d = shiftU(D) - D; % forward
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
  
%-- whole matrix derivatives
function shift = shiftD(M)
  shift = shiftR(M')';

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
  shift = shiftL(M')';
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    

  





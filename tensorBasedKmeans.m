function segmentedImage = tensorBasedKmeans(I,k,form)

  I = double(rgb2gray(I));  
  g = fspecial('Gaussian');
  I = imfilter(I, g);
  % computing the structure tensors
  % getting Ix and Iy
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
    sigma = 0.8;
    hsize = 3;
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
     
%COMPUTATION OF STRUCTURE TENSORS DONE

numberOfFeatures = 4; %as every tensor represented by 4 values
% centroid init randomly
% each row of centroidMatrix represents one centroid
centroidMatrix = rand(k,numberOfFeatures);
for i=1:k
    min1 = min(min(tensorI1));
    centroidMatrix(i,1) = centroidMatrix(i,1)*(max(max(tensorI1))- min1) + min1;
     min2 = min(min(tensorI2));
    centroidMatrix(i,2) = centroidMatrix(i,2)*(max(max(tensorI2))- min2) + min2;
     min3 = min(min(tensorI3));
    centroidMatrix(i,3) = centroidMatrix(i,3)*(max(max(tensorI3))- min3) + min3;
     min4 = min(min(tensorI4));
    centroidMatrix(i,4) = centroidMatrix(i,4)*(max(max(tensorI4))- min4) + min4;
end
%CENTROID RANDOM INITIALISATION DONE

%assigning centroids to points
%taking out every point's distance from the centroids

clusterAssignments = zeros(size(I,1),size(I,2));
for inc = 1: 5
  
  % finalClusterAssignments pixels to closest centroid
  for i = 1 : size(I,1)
      for j=1:size(I,2)
    %initialising the minimum distance of this (i,j)th point from a centroid
    %to it's distance from the first centroid in the centroidMatrix
        minDistanceToACentroid = distance(reshape(STMatrix(i,j,:,:),2,2),reshape(centroidMatrix(1,:),2,2),form);
        clstr = 1;
        for c = 2 : k
         distanceToACentroid = distance(reshape(STMatrix(i,j,:,:),2,2),reshape(centroidMatrix(c,:),2,2),form);
          if( minDistanceToACentroid >= distanceToACentroid)
            clstr = c;
            minDistanceToACentroid = distanceToACentroid;
          end
        end

    %book-keeping for what cluster has been finalClusterAssignmentsed to each pixel
    clusterAssignments(i,j) = clstr;
      end
  end

  
  % updating centroids by recalculating their positions

    for c=1:k
        newC = meanFunc(STMatrix,clusterAssignments,c,form);
        centroidMatrix(c,:) = reshape(newC,1,4);
    end
  end
  %updating the centroids by the new means
  
  segmentedImage = clusterAssignments;
  imagesc(segmentedImage);
  
%   segment = zeros(size(I,1),size(I,2));
%   
%   for i=1:k
%     figure(i);
%         [x,y] = find(segmentedImage == i);
%         for l=1:length(x)
%             segment(x(l),y(l))=1;
%         end
%     imagesc(segment);
%   end
%   for i=1:k
%       [x,y] = find(segmentedImage == i);
%       for l=1:length(x)
%        finalSegments(l,x,y)= 1;
%       end
%   end
  
  
end


%computing means based on the distance type
function u = meanFunc(st,clusterAssignments,clusterValue,form)
  [ptsx,ptsy] = find(clusterAssignments == clusterValue);
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
end  
  
function dist = distance(st,u,form)
 dist = 0;
   %distance of a point from a centroid rep by u
   if form == 1    
        for i=1:2
            for j=1:2
            dist = dist + (st(i,j)-u(i,j))^2;    
            end
        end
   elseif form == 2
          dist = getJDist(st,u);
   end
end

function d2 = getJDist(tensor,u)
      d2 = 1/2*trace(inv(tensor)*u +inv(u)*tensor)-2;
end
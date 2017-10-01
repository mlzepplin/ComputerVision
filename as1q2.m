clear
imageFile = '\\Client\C$\Users\rishabh\Desktop\shelf2.jpeg';

%getting 2 vanishing points to get the equation of vanishing line
v1 = vanishingPoint(imageFile);
v2 = vanishingPoint(imageFile);

%vanishing line is the line that passes through both v1 and v2
%polyfit([x1,x2],[y1,y2])gets corresponding line (y=ax+b), a = coefficients (1),b1 = coefficients (2)
coefficients = polyfit([v1(1,1),v2(1,1)],[v1(1,2),v2(1,2)], 1);

%as the euclidean line cx + dy + e =0 represents(is a projection of) the plane with normal
%(c,d,e) of homogenous coordinates
l = [coefficients(1),-1,coefficients(2)]; %vanishing line

%getting the vertical vanishing point of (t1,b1) and (t2,b2)
%INPUT AS [t1,b1,t2,b2]
points = Get2DPoints(imageFile, 4);

t1 = points(:,1);
b1 = points(:,2);

t2 = points(:,3);
b2 = points(:,4);

v = vanishingPointForGivenRays([t1',1],[b1',1],[t2',1],[b2',1]);


%getting vanishing point corresponding to b1 and b2
u = cross(cross([b1',1],[b2',1]),l);


%computing t~1 i.e. t tilde 1, here a.k.a tp1
l2 = cross(v,[b2',1]);
tp1 = cross(cross([t1',1],u),l2);




%getting equation of l1
%polyfit([x1,x2],[y1,y2])gets corresponding line (y=ax+b), a = coefficients (1),b1 = coefficients (2)
coefficients = polyfit([t1(1,1),b1(1,1)],[t1(2,1),b1(2,1)], 1);
%slolpe of l1 be m
m = coefficients(1);

%l3 is parallel to l1 and passing throught b2
%y-y1 = m(x-x1) => mx -y + (y1-mx1) = 0
l3 = [m,-1,b2(2,1)-m*b2(1,1)];

%obtaining distances i.e. the ^ values
tp1Cap = cross(cross([t1',1],tp1),l3);
t2Cap = cross(cross([t1',1],[t2',1]),l3);

g1 = v - t2Cap;
g2 = v - tp1Cap;

distRatio = abs(dot(tp1Cap,g1)/dot(t2Cap,g2))

 
function v = vanishingPoint(imageFile)

%line p1 made by euclidean ponits p11 and p12
p1Points = Get2DPoints(imageFile, 2)

%(parallel in the real world) line p2 made euclidean points p21 and p22
p2Points = Get2DPoints(imageFile, 2)

%homogenising and converting into rays, and into  row vectors(to comply with matlab cross func)
p1PointsPrime = vertcat(p1Points, ones(1,2));
p2PointsPrime = vertcat(p2Points, ones(1,2));

r1 = [p1PointsPrime(:,1)']
r2 = [p1PointsPrime(:,2)']

r3 = [p2PointsPrime(:,1)'];
r4 = [p2PointsPrime(:,2)'];

vPrime = vanishingPointForGivenRays(r1,r2,r3,r4);

%getting the euclidean vanishing point v, note the last column of the
%vanishing point will be zero in homogenoue(Prime)representation
v = vPrime(1,1:2);
    
return ;
end

function vPrime = vanishingPointForGivenRays(r1,r2,r3,r4)
%returns the homogeneuos rep of the vanishing point
%computing the normal that defines the plane formed by the 2 rays r1 and r2
n1 = cross(r1,r2);
%similarly for n2
n2 = cross(r3,r4);

%intersection of these two planes will give one ray
%and when that ray's point of intersection on the image plane 
%will be the vanishing point for these set of parallel lines
vPrime = cross(n1,n2);



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

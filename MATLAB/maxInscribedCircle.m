function [results, center,radius]  = maxInscribedCircle(data1)
%Prem Rachakonda, 2018

xx = data1(:,1);
yy = data1(:,2);
%[vx,vy] = voronoi(xx,yy); 
vxy1 = voronoin(data1); %voronoin is multidimensional and outputs unique rows

%The voronoi points are unique, but the output of voronoi() are the ends 
% of the edges, so unique points (nodes) have to be extracted. 
%vxy1 = unique([vx(:) vy(:)],'rows');
vx1 = vxy1(:,1);
vy1 = vxy1(:,2);

%Find the points of the voronoi points that are within a convex hull. The
%MIC has to be within this convex hull, so the "search" for MIC doesn't
%need to extend beyond this envelope
k = convhull(xx,yy);
idx1 =inpolygon(vx1,vy1, xx(k),yy(k));

%These are the voronoi nodes within the envelope of the circle
vx2 = vx1(idx1);
vy2 = vy1(idx1);


%For each voronoi point, find the distances to all the points
%Find the radial distances of vornoi points (vx2,vy2) from the input points
%Then find the minimal radius at each node
for ii = 1:length(vx2)
    radiusAtNode= (sqrt( bsxfun(@minus,xx,vx2(ii)).^2 + bsxfun(@minus,yy,vy2(ii)).^2));
    radMin(ii) = min(radiusAtNode);
end
%Each row represents one voronoi node

%For all the circles at all the nodes, find the largest circle radius
[radius,idxa] = max(radMin);
center = [vx2(idxa),vy2(idxa)];

results.Center = center;
results.Radius = radius; 
function resids1 = findDirResids(fitCenter,trueRadius,data1)
% fitCenter should by a 1x3 row vector
% trueRadius should be a single scalar value
% data1 should be a Nx3 matrix (N rows of 3D point coordinates)
%
% Refer to "Fitting Spheres to Range Data From 3-D Imaging Systems" by Marek Franaszek, Geraldine S. Cheok, Kamel S. Saidi, and Christoph Witzgall
% http://www.nist.gov/customcf/get_pdf.cfm?pub_id=861541
%
% Author: Prem Rachakonda
% Date: 2/22/2015
% Version: 1.0

%%Sanitize the data
[m1,n1] = size(fitCenter)
if m1==3 & n1 ==1
    fitCenter = fitCenter';
else
    error('fitCenter should be a row vector of 1x3 size');
end

[m,n] = size(data1);
if(n<3 | m<4)
    data = data';
end
[m,n] = size(data1);
if(n<3 | m<4)
    error('Not enough data');
end

%% Formulate a matrix of centers to make it easy for division later on
centerMat = repmat(fitCenter,length(data1),1);
RR = trueRadius;

%% Find the radial distances of each points and formulate a matrix to make it easy for division later on.
rj = rssq2(data1,2);
rjMat = repmat(rj,1,3);

%%Generate unit vectors for each data point.
uData1 = data1./rjMat;

%% According to Markek F's paper, find the dot product (pj) and the magnitude of the cross product(qj)
pj = dot(uData1,centerMat,2);
qj = rssq2( cross(uData1,centerMat) ,2);

%% Based on whether qj<RR or qj>=RR, find the residuals
idx = qj<RR;
res1 = pj-rj - sqrt(RR^2 -qj.^2); % Vector corresponding to the point intersects the sphere
res2 = sqrt((pj-rj).^2 + (qj-RR).^2); % Vector corresponding to the point does not intersect the sphere

resids1(idx) = res1(idx);
resids1(~idx) = res2(~idx);
resids1 = resids1'; %For some reason, if res1 is column vector, resids1 is a row vector?? 

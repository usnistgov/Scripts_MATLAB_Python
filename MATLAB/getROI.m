function data2 = getROI(data1,center1,dist1)
%Truncate datasets based on a center and the distance of points from the center
minVal = center1-dist1;
maxVal = center1+dist1;
idx1 = data1(:,1)>minVal(1) & data1(:,1)<maxVal(1);
idx2 = data1(:,2)>minVal(2) & data1(:,1)<maxVal(2);
idx3 = data1(:,3)>minVal(3) & data1(:,1)<maxVal(3);
idxA = idx1&idx2&idx3;
data2 = data1(idxA,:);
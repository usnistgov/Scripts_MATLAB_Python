function [dataFinal,rng2] = firstNPoints(data1,N)
%Calculates first N points on a sphere from the center
[~,~,rng1] = cart2sph(data1(:,1),data1(:,2),data1(:,3));
[rng2,idx] = sort(rng1);
data2 = data1(idx,:);
dataFinal = data2(1:N,:);

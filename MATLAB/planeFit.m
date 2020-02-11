function result = planeFit3(data)
centroid = mean(data,1);
dataC = bsxfun(@minus,data,centroid);
[U,S,V] = svd(dataC,0);
P = V(:,3);
P = P/sqrt(sum(P.*P));
d = -dot(P,cm);
P(4) = d;

P = P*sign(P(1)); %Just making sure that the sign of the equations are always the same
a = P(1); b = P(2); c= P(3); d = P(4);
result.Parameters = [a b c d];

function result = planeFit2(data)
% result = planeFitT(data)
% data has to be of the form N rows, 3 columns for this to work. 
% 
% Prem Rachakonda 

%Translate the data to its centroid
centroid = mean(data,1);
dataC = bsxfun(@minus,data,centroid);

%Calculate the eigen vectors and obtain the normal vector
[V,~] = eig(dataC'*dataC);
P = V(:,1);

P = P/sqrt(sum(P.*P));
d = -dot(P,centroid);
P(4) = d;

P = P*sign(P(1)); %Just making sure that the sign of the equations are always the same
a = P(1); b = P(2); c= P(3); d = P(4);
xx= data(:,1); yy = data(:,2); zz = data(:,3);

result.Parameters = [a b c d];
result.Residuals = (a*xx + b*yy + c*zz + d)/sqrt(a^2 + b^2 + c^2);
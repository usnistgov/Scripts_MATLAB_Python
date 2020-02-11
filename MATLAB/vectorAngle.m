function vAngle = vectorAngle(vec1, vec2)
% vAngle = vectorAngle(vec1, vec2)
% 
% Author:  Prem Rachakonda
% Date:    April 02, 2016 
% 
% Function Description: 
%         This function takes two vectors and calculates the angle between
%         them. 
% Input Parameters: 
%         1. vec1:  First vector (Nx3 array)
%         2. vec2: Second vector (Nx3 array)
%   
% Returned Values:  
%         1. vAngle: Angle between the two input vectors in radians (Nx1)
%

sinValue = rssq2n(cross(vec1,vec2),2);
cosValue = dot(vec1',vec2')';
vAngle = atan2(sinValue,cosValue);

function result = rssq2n(data1,DIM)
result = data1.^2;
result = sqrt(sum(result,DIM));

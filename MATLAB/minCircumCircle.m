function [results, center1,radius1] = minCircumCircle(data1)
%Prem Rachakonda, 2018

%We first setup the objective function. The initialGuess is the variable to
%be minimized and its initial value is the centroid.
objectiveFunction = @(initGuess)maxRadialDist(initGuess,data1);
[center1,radius1] = fminsearch(objectiveFunction,mean(data1),data1);
%x0 = center1 
%fval = radius1 (function value at the minmized value

results.Center = center1;
results.Radius = radius1; 

function dist1  = maxRadialDist(initGuess,data1)
%This is the cost function that fminsearch() has to minimize
dist1 = max( rssq2( bsxfun(@minus,data1,initGuess), 2) );



function result = sphereFitLSQ1(data,dummyVal1)
% This function takes 3 dimensional data of a sphere in a Cartesian 
% coordinate system and returns the paramters of a sphere (Center and
% radius)
%
%The input data is in the form of Nx3 matrix
%The output is a structure with the sphere parameters.
%The function also takes in an optional parameter "dummyVal1" which does
%not do anything, and is retained for consistency purposes with other
%similar algorithms

%If no "dummyVal" is given, just assign some random number. It is not used anyway. 
if nargin<2
    dummyVal = 1;
end


%First get an initial estimate based on a linearized form of the sphere equation
[ro,ao,bo,residualo] = sphereFitInitGuess(data);
initGuess = [ao bo ro];

%Now formulate the objective function and the options for lsqnonlin() 
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunEvals',1500,'TolFun',1E-9,'display','off');
objectiveFunction = @(initGuess)calcResiduals(initGuess,data);
[result3] = lsqnonlin(objectiveFunction,initGuess,[],[],options);

%The above part works fine. The part below tries to use "multiple start" to
%make sure that the LMA algo didn't get stuck at the local minimum
if exist('createOptimProblem')>0
    options2 = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
    problem = createOptimProblem('lsqnonlin','x0',initGuess,'objective',objectiveFunction, 'options', options2);
    %ms = MultiStart('PlotFcns',@gsplotbestf);
    ms = MultiStart('Display','off');
    [result3,f] = run(ms,problem,20);
end

%We have other algorithms where the output parameters have 4 elements, so
%to keep it consistent, the 4th parameter is added.
if(length(result3)<3)
    result3(3) = trueRadius;
end

% The "residuals" and "Residuals" are the same parameters, but retained for
% legacy purposes
result.Parameters = result3;
result.residuals = calcResiduals(result3,data);
result.Residuals = result.residuals;
result.Center = result3(1:end-1)';
result.Radius = result3(end);



function residuals = calcResiduals(initGuess,data)
%Function to calculate the residuals based on the sphere parameters and the
%data.

ao = initGuess(1);
bo = initGuess(2);
ro = initGuess(end);

xx = data(:,1);
yy = data(:,2);
residuals = sqrt( (xx-ao).^2 + (yy-bo).^2) - ro;



function [r,a,b,residual] = sphereFitInitGuess(data)
%This is a function to calculate the parameters of a sphere as an initial
%estimate. 
xx = data(:,1);
yy = data(:,2);


AA = [-2*xx, -2*yy , ones(size(xx))];
BB = [ -(xx.^2+yy.^2)];

YY = mldivide(AA,BB);

a = YY(1);
b = YY(2);
D = YY(3); % a^2 + b^2 + c^2 -r^2(where a,b,c are centers)
r = sqrt((a^2+b^2)-D);
residual = AA*YY-BB;


% [x2,y2,z2]=sphere(128);
% x2=r*x2+a; %Multiply the x,y,zs with 'r' to scale the unit circle and offset it by adding a,b,c (the center)
% y2=r*y2+b;
% z2=r*z2+c;
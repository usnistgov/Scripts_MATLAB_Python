function X = intersectPlanes(params1, params2, params3, params4)
%Intersect 4 planes and return the intersection point
A = [params1(1:3); params2(1:3); params3(1:3); params4(1:3)];
B = [-params1(4); -params2(4); -params3(4); -params4(4)];
X = A\B; X = X';

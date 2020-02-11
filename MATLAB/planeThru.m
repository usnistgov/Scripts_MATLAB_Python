function planeParams = planeThru(p1, p2, p3)
%Calculate plane equation through 3 points
normal = cross(p1 - p2, p1 - p3);
d = p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3);
normal(4) = -d;
planeParams = normal;

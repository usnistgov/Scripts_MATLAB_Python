function data1 = createRandomSphere(a,b,c,r)

if nargin <1
    a = 1;
    b = 2;
    c = 3;
    r = 10;
end

[x,y,z] = sphere(100);
x1 = a + r*x(:);
y1 = b + r*y(:);
z1 = c + r*z(:);

x2 = x1 + randn(size(x1))*r/100;
y2 = y1 + randn(size(y1))*r/100;
z2 = z1 + randn(size(z1))*r/100;
plot3(x2,y2,z2);
data1 = [x2 y2 z2];

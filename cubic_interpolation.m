%cubic interpolation
%三次插值

% out = 4* epsi* ((sigma/ x)^12 - (sigma/ x)^6 );

epsi = 1;
sigma = 1;

a = 0.98810055;
b = 1.74615763;

f1 = LJ(epsi,sigma,a);
g1 = grad(epsi,sigma,a);

f2 = LJ(epsi,sigma,b);
g2 = grad(epsi,sigma,b);

c1 = (g1+ g2- 2* (f2- f1)/ (b- a))/(b-a)^2;
c2 = 3* (f2- f1)/(b- a)^2 - (2* g1+ g2)/(b- a);
c3 = g1;
c4 = f1;

extention = 0.1;

%画图
x = a-extention:(b-a)/50:b+extention;
% x = a:(b-a)/50:b;
y1 = LJ(epsi,sigma,x);

y2 = cubicInterpolation(x,a,c1,c2,c3,c4);




plot(x, y1, x, y2, 'LineWidth', 1);
xlabel('x2');
grid on;




function out = LJ(epsi,sigma,x)
    out = 4* epsi* ((sigma./ x).^12 - (sigma./ x).^6 );
end

function out = grad(epsi,sigma,x)
    out = 4* epsi* (-12* sigma^12./ x.^13 + 6* sigma^6./ x.^7 );
end

function out = cubicInterpolation(x,x0,c1,c2,c3,c4)
    out = c1* (x-x0).^3 + c2* (x-x0).^2 + c3* (x-x0) + c4;
end





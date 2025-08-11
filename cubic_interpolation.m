%cubic interpolation
%三次插值

clear %清理变量
clf %清理图像

% out = 4* epsi* ((sigma/ x)^12 - (sigma/ x)^6 );

%力场参数
epsi = 1;
sigma = 1;

%三阶插值近似求解的步长和精度
% a = 0;
% b = 0;
change = 0;%函数值交换用的中间值
accu = 0.01;
s = 0;
z = 0;
w = 0;

t1 = 1.13804235;
t2 = 0.74996936;

f1 = LJ(epsi,sigma,t1);
g1 = grad(epsi,sigma,t1);

f2 = LJ(epsi,sigma,t2);
g2 = grad(epsi,sigma,t2);

c1 = (g1+ g2- 2* (f2- f1)/ (t2- t1))/(t2-t1)^2;
c2 = 3* (f2- f1)/(t2- t1)^2 - (2* g1+ g2)/(t2- t1);
c3 = g1;
c4 = f1;


extention = 0.1;

x = t1+ sign(g1)* extention:(t2-t1)/50:t2- sign(g1)* extention;
y1 = LJ(epsi,sigma,x);

y2 = cubicInterpolation(x,t1,c1,c2,c3,c4);

plot(x, y1, x, y2, 'LineWidth', 1);
xlabel('x2');
grid on;



condition_t = t2 < t1;

s = 3* (f2 - f1)/ (t2 - t1);
z = s - g1 - g2;
w = sqrt(z*z - g1*g2);

tNew = (1-condition_t)*t1 + condition_t*t2 +...
    (-1)^condition_t * (t2-t1)*...
    (1- ( (1-condition_t)*g2 + condition_t*g1 + w+ z)/...
    ((-1)^condition_t*(g2- g1) + 2*w));

fNew = LJ(epsi,sigma,tNew);
gNew = grad(epsi,sigma,tNew);
fprintf('求解出的新点：a1: %.8f, f1: %.4f, g1: %.4f, \n',tNew,fNew,gNew);
hold on
plot(tNew,cubicInterpolation(tNew,t1,c1,c2,c3,c4),'kx','MarkerSize',12)


function out = LJ(epsi,sigma,x)
    out = 4* epsi* ((sigma./ x).^12 - (sigma./ x).^6 );
end

function out = grad(epsi,sigma,x)
    out = 4* epsi* (-12* sigma^12./ x.^13 + 6* sigma^6./ x.^7 );
end

function out = cubicInterpolation(x,x0,c1,c2,c3,c4)
    out = c1* (x-x0).^3 + c2* (x-x0).^2 + c3* (x-x0) + c4;
end



% 报废代码
% 
% condition_t = 0;
% s = 3* (f2 - f1)/ (t2 - t1);
% z = s - g1 - g2;
% w = sqrt(z*z - g1*g2);
% tNew = t1 + (t2-t1)* (1- (g2+ w+ z)/(g2- g1+ 2*w));



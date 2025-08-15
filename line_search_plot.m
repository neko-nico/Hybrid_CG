%visualization of line search
%固定搜索方向后，画一下沿着搜索方向的线搜索迭代图

clear %清理变量
clf %清理图像

%画图参数
set(groot , 'DefaultAxesFontName', '宋体', 'DefaultAxesFontSize', 15);  % 设置坐标轴标签的默认字体和大小
set(groot, 'defaultTextFontName', 'Times New Roman', 'defaultTextFontSize', 14);
set(groot, 'defaultLineMarkerSize', 10);

%set(0, 'DefaultLineLineWidth', 2);  % 设置默认线条的粗细
%set(groot , 'DefaultLineMarkerSize', 12);  % 设置默认点的半径

epsi = 1;
sigma = 1;

%极小值点
min_x = 2^(1/6)* sigma;

%三阶插值近似求解的步长和精度
a = 0;
b = 0;
condition_t = 0;
change = 0;%函数值交换用的中间值
accu = 0.1;
s = 0;
z = 0;
w = 0;

%三阶插值图像绘制
c1 = 0;
c2 = 0;
c3 = 0;
c4 = 0;

%c初始点

t1 = 0.92 + (1.2 - 0.92) * rand;
% t1 = 0.96;
% fprintf('t1 = %.8f;',t1)
f1 = LJ(epsi,sigma,t1);
g1 = grad(epsi,sigma,t1);
% fprintf('t1: %.8f, f1: %.4f, g1: %.4f, \n',t1,f1,g1);

%步长
h = abs(g1)/100;
% h = 2*f1/g1;

%循环线限制
limit = 0;
limit_max = 2;

%subplot
extention = 0.001;
sub_cols = 2;
sub_rows = ceil(limit_max/sub_cols);

while (abs(g1) > accu && limit < limit_max)|| limit == 0
    limit = limit +1;
    fprintf('\n第 %d 个点,\n',limit);
    fprintf('t1 = %.8f;\n',t1)
    fprintf('h = %.8f;\n',h)
    fprintf('t1: %.8f, f1: %.4f, g1: %.4f, \n',t1,f1,g1);
    
    if g1 < 0
        h = abs(h);                                                                           
    else
        h = - abs(h);
    end

    t2 = t1 + h;
    f2 = LJ(epsi,sigma,t2);
    g2 = grad(epsi,sigma,t2);
    fprintf('t2: %.8f, f2: %.4f, g2: %.4f,\n',t2,f2,g2)

    if abs(g2) < accu %&& abs(g2) < abs(g1)
        t1 = t2*2^sign(g1*g2);
        f1 = f2;
        g1 = g2;
        break;
    end
    
    if g1 * g2 < 0

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
        fprintf('求解出的新点：tNew: %.8f, fNew: %.4f, gNew: %.4f, \n',tNew,fNew,gNew);

        %绘制图像
        c1 = (g1+ g2- 2* (f2- f1)/ (t2- t1))/(t2-t1)^2;
        c2 = 3* (f2- f1)/(t2- t1)^2 - (2* g1+ g2)/(t2- t1);
        c3 = g1;
        c4 = f1;

        x = t1+ sign(g1)* extention:(t2-t1)/1000:t2- sign(g1)* extention;
        y1 = LJ(epsi,sigma,x);
        y2 = cubicInterpolation(x,t1,c1,c2,c3,c4);

        subplot(sub_rows,sub_cols,limit)
        plot(x, y1, x, y2, 'LineWidth', 1);
        % xlabel('x');
        % ylabel('y');
        title(['第 ', num2str(limit), ' 个点'])
        grid on;

        %绘制箭头
        hold on;
        plot([t1, t2], [f1, f2], 'g-', 'LineWidth', 1.5);
        text((t1+t2)/2, (f1+f2)/2,['h = ',num2str(h)]);
        
        %quiver(t1,f1,tNew-t1,fNew-f1, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1/norm([tNew-t1,fNew-f1]));
        plot([t1, tNew], [f1, fNew], 'm-', 'LineWidth', 1.5);
        text((t1+tNew)/2, (f1+fNew)/2,['\Deltat = ',num2str(tNew-t1)]);
        fprintf('当前迭代 Delta t = %.8f\n',tNew-t1);
        % legend('LJ势能曲线','cubic拟合曲线','搜索点迭代')
        % hold off;

        if ( t1<tNew && tNew<t2 )||( t1>tNew && tNew>t2 )
            %不允许跨过极小值点
            if gNew * g1 > 0 || abs(gNew) < abs(g1)/2
                t1 = tNew;
                f1 = fNew;
                g1 = gNew;
                fprintf('取了新的点！\n')
            else
                fprintf('没有取新的点…\n')
                h = h*1.5;%少收缩一些
            end
                
            if abs(gNew) > accu || limit == 1
                h = h/6;
                fprintf('步长 h 太长了！！ %.5f\n',h);
           
                h = -sign(g1)*abs(h);
                plot([t1, t1+h], [f1, LJ(epsi,sigma,t1+h)], 'r-', 'LineWidth', 1.5);
                text(t1+h/2, (f1+LJ(epsi,sigma,t1+h))/2,['h = ',num2str(h)]);
                fprintf('下一轮迭代 h = %.8f\n',h);
                plot(min_x,LJ(epsi,sigma,min_x),'kx')
                legend('LJ势能曲线','cubic拟合曲线','初始区间','当前迭代','下一轮迭代')
                hold off

            else

                legend('LJ势能曲线','cubic拟合曲线','初始区间','当前迭代')
            end


        else
            h = h/4;
            t1 = (t1+t2)/2;
            f1 = LJ(epsi,sigma,t1);
            g1 = grad(epsi,sigma,t1);
            fprintf('三次插值出错，取中点\n')
        end

        
    elseif abs(g2) < 10*abs(g1) %|| f2 < f1
        
        %画图
        x = t1-extention:3*h/50:t1+3*h+extention;
        y1 = LJ(epsi,sigma,x);
        
        subplot(sub_rows,sub_cols,limit)
        plot(x, y1, 'LineWidth', 1);
        % xlabel('x');
        % ylabel('y');
        title(['第 ', num2str(limit), ' 个点'])
        grid on;
        
        %绘制箭头
        hold on;
        plot([t1, t2], [f1, f2], 'g-', 'LineWidth', 1.5);
        text((t1+t2)/2, (f1+f2)/2,['h = ',num2str(h)]);

        %quiver(t1,f1,tNew-t1,fNew-f1, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1/norm([tNew-t1,fNew-f1]));
        plot([t1, t2], [f1, f2], 'm-', 'LineWidth', 1.5);
        text((t1+t2)/2, (f1+f2)/2,['h = ',num2str(h)]);
        fprintf('当前迭代 h = %.8f\n',h);

        h = 2*h;
        t1 = t2;
        f1 = f2;
        g1 = g2;
        fprintf('前进,而且步长 h 太短了…… %.5f\n',h);

        plot([t1, t1+h], [f1, LJ(epsi,sigma,t1+h)], 'r-', 'LineWidth', 1.5);
        text(t1+ h/2, (f1+ LJ(epsi,sigma,t1+h))/2,['h = ',num2str(h)]);
        fprintf('下一轮迭代 h = %.8f\n',h);
        plot(min_x,LJ(epsi,sigma,min_x),'kx')

        legend('LJ势能曲线','初始区间','当前迭代','下一轮迭代')
        hold off
    else
        h = h/6;
        fprintf('梯度异常，返回\n');
    end

end


function out = LJ(epsi,sigma,x)
    out = 4* epsi* ((sigma./ x).^12 - (sigma./ x).^6 );
end

function out = grad(epsi,sigma,x)
    out = 4* epsi* (-12* sigma^12./ x.^13 + 6* sigma^6./ x.^7 );
end

function out = cubicInterpolation(x,x0,c1,c2,c3,c4)
    out = c1* (x-x0).^3 + c2* (x-x0).^2 + c3* (x-x0) + c4;
end



%绘制L-J势能图像 - ok
%
% epsi = 1;
% sigma = 1;
% x1 = 0;
% 
% % 生成 x1 的值
% x2 = linspace(1, 2, 100);  % 避免除以 0，所以不要包含 x1 = x2
% 
% % 计算 f 对应的值
% fval = arrayfun(@(x) LJ(epsi, sigma, x, x1), x2);
% 
% % 绘图
% plot(x2, fval, 'LineWidth', 2);
% xlabel('x2');
% ylabel('f(x2)');
% title('Lennard-Jones Potential');
% grid on;
% 
% function out = LJ(epsi,sigma,r1,r2)
%     sigmarj = sigma/(r1-r2);
%     out = 4* epsi* (sigmarj^12-sigmarj^6);
% end
% 
% function out = LJ(epsi,sigma,x)
%     out = 4* epsi* ((sigma/ x)^12 - (sigma/ x)^6 );
% end
% 
% function out = grad(epsi,sigma,x)
%     out = 4* epsi* (-12* sigma^12/ x^13 + 6* sigma^6/ x^7 );
% end

%嵌入原子势

clear
set(0, 'DefaultAxesFontName', '宋体', 'DefaultAxesFontSize', 18);  % 设置坐标轴标签的默认字体和大小
set(0, 'DefaultTextFontSize', 18);  % 设置文本的默认字体大小

%体系维数+粒子数目
dimension = 3;
pointsNum = 100;
n = dimension*pointsNum;%总维数

%力场参数
epsilon = 1;
sigma = 3.46;

%初始粒子坐标取值范围
iniRange = sigma* 2^(1/6)* (0.5+ (3*pointsNum/ (4* pi* sqrt(2)))^ (1/3) );

% 生成随机点的空坐标
pointsList = zeros(n,1);

%第一个点要在圆内
randomPoint = (2*rand(dimension,1)-1)*iniRange;
while norm(randomPoint) > iniRange
    randomPoint = (2*rand(dimension,1)-1)*iniRange;
end
pointsList(1:dimension) = randomPoint;

for i = 2:pointsNum
    randomPoint = (2*rand(dimension,1)-1)*iniRange;
    minDis = distence(pointsList,randomPoint,dimension);
    while norm(randomPoint) > iniRange || minDis(1) < sigma
        randomPoint = (2*rand(dimension,1)-1)*iniRange;
        minDis = distence(pointsList,randomPoint,dimension);
        % fprintf(logfid, '太近了，重开\n');
    end
    pointsList((i-1)*dimension+1:i*dimension) = randomPoint;
end

fprintf('initial\n')

load('E:\MatlabProjects\hybrid_CG\Result\0909-164514_pointsList_1_check.mat');

energy = potential(epsilon,sigma,pointsList,dimension,pointsNum);
gr = gradient(epsilon,sigma,pointsList,dimension,pointsNum);%初始点梯度
dr = - gr;%初始搜索方向为负梯度方向

load('E:\MatlabProjects\hybrid_CG\Result\0909-164514_dr_1_check.mat');

%取样点
sample = 100;

%用数组记录grad,direction,energy,h,方向变化间隔
k = 1;
t = 1;
glist = gr;
dlist = dr;
elist = energy;
%搜索步长
% h = max(abs(gr))/(100*sample);
h = 0.15447/sample;
%fprintf('初始步长: %.4f\n',h)
% hlist = h;
% section = 1;
nantest = 0; %检测NaN报错
maxMove = 0.3;

drgrList = dr'* gr;
plist_hdr = h*dr;


for i = 0:sample
    pointsList = pointsList + plist_hdr;

    energy = potential(epsilon,sigma,pointsList,dimension,pointsNum);
    gr = gradient(epsilon,sigma,pointsList,dimension,pointsNum);

    elist = [elist,energy];
    drgrList = [drgrList, dr'* gr];
end

x = 1:length(elist);

f = figure('Position', [740, 680, 800, 600]);
% 绘制
yyaxis left
plot(x, elist, '-', 'LineWidth', 1.5, 'MarkerSize', 5);
ylabel('函数值');

yyaxis right
plot(x, drgrList, '-', 'LineWidth', 1.5, 'MarkerSize', 5);
ylabel('方向导数');
yline(0, '--k', 'LineWidth', 1);  % '--k' 表示黑色虚线

xlabel('取样点');
legend('函数值','方向导数')


title('EAM一维势能面');


ax = gca;
% 改右轴字体颜色
yyaxis right
ax.YColor = [0 0 0];
% 改左轴字体颜色
yyaxis left
ax.YColor = [0 0 0];
% 改横坐标字体颜色
ax.XColor = [0 0 0];
grid on;




%粒子群中，离某个位置最近的的原子 的序号和距离
function out = distence(pointsList,vector,dim)
    pointi = 1;
    min = norm(vector - pointsList(1:dim));
    if min == 0
        pointi = 2;
        min = norm(vector - pointsList(dim+1:2*dim));
    end

    for i = dim+1:dim:length(pointsList)
        if min > norm(vector-pointsList(i:i+dim-1)) && norm(vector-pointsList(i:i+dim-1)) ~= 0
            min = norm(vector-pointsList(i:i+dim-1));
            pointi = (i+dim-1)/dim;
        end
    end
    out = [min,pointi];
end

%记录所有粒子间距离的数组
function out = averageDistence2List(pointsList,dim,num)
    out = zeros(num,1);
    for i = 1:num
        mid = 0;
        ri = pointsList((i-1)*dim+1:i*dim);
        for j = 1:num
            mid = mid + norm(ri - pointsList((j-1)*dim+1:j*dim))^2;
        end
        out(i) = mid/num;
    end
end

%势能
function out = potential(epsi,sigma,plist,dim,num)
out = 0;
for i = 0:num-2
    ri = plist(i*dim+1:(i+1)*dim);
    for j = i+1:num-1
        sigmarj = sigma/norm(ri-plist(j*dim+1:(j+1)*dim));
        out = out + 4* epsi* (sigmarj^12-sigmarj^6);
    end
end
end



%梯度
function out = gradient(epsi,sigma,plist,dim,num)
out = zeros(dim*num,1);
mid0 = zeros(dim,1);
for i = 0:num-1
    ri = plist(i*dim+1:(i+1)*dim);
    mid = mid0;
    for j = 0:num-1
        if j ~= i
            rirj = ri - plist(j*dim+1:(j+1)*dim);
            mid = mid + 4* epsi* rirj * (-12*sigma^12/norm(rirj)^14 + 6*sigma^6/norm(rirj)^8);
        end
    end
    out(i*dim+1:(i+1)*dim) = mid;
end
end


%将粒子坐标输出成Vesta文件格式
function writeVestaFile(filename,pointsList,dimension,pointsNum)

    %写入文件需要的杂七杂八的内容
    text0 = fileread("origin0.txt");
    text1 = fileread("origin1.txt");
    text2 = fileread("origin2.txt");
    text3 = fileread("origin3.txt");
    text4 = fileread("origin4.txt");
    fileID = fopen(filename,'w');
    
    fprintf(fileID,text0);
    cellP = 1:dimension:dimension*pointsNum;
    cellP = [max(pointsList(cellP))-min(pointsList(cellP)),...
             max(pointsList(cellP+1))-min(pointsList(cellP+1)),...
             max(pointsList(cellP+2))-min(pointsList(cellP+2)),];
    fprintf(fileID,' %9.6f  %9.6f  %9.6f  ',cellP);
    fprintf(fileID,text1);
    
    pointsList = 0.98*(pointsList-min(pointsList))/(max(pointsList)-min(pointsList));
    
    for i = 1:pointsNum
        fprintf(fileID,'  %d  C         C%d  1.0000   %1.6f   %1.6f   %1.6f    1a       1\n',i,i,pointsList(dimension*i-2),pointsList(dimension*i-1),pointsList(dimension*i));
        fprintf(fileID,'                            0.000000   0.000000   0.000000  0.00\n');
    end
    
    fprintf(fileID,'  0 0 0 0 0 0 0\n');
    fprintf(fileID,'THERI 1\n');
    
    for i = 1:pointsNum
        fprintf(fileID,'  %d         C%d  0.050000\n',i,i);
    end
    
    fprintf(fileID,text2);
    
    for i = 1:pointsNum
        fprintf(fileID,'  %d         C%d  %5.4f 128  73  41 128  73  41 204  0\n',i,i,3/pointsNum^(1/dimension));
    end
    
    fprintf(fileID,text3);
    fprintf(fileID,'%1.4f',5/pointsNum^(1/dimension));
    fprintf(fileID,text4);
    
    fclose(fileID);
end



%嵌入原子势

clear
set(0, 'DefaultAxesFontName', '宋体', 'DefaultAxesFontSize', 18);  % 设置坐标轴标签的默认字体和大小
set(0, 'DefaultTextFontSize', 18);  % 设置文本的默认字体大小

%体系维数+粒子数目
dimension = 3;
pointsNum = 100;
n = dimension*pointsNum;%总维数

%力场参数
r0Prm = 2.556;
rcPrm = 5;
E0Prm = 3.54;
Phi0Prm = 4.86;
alphaPrm = 5.9;
betaPrm = 5.85;
gammaPrm = 6.3;

%初始粒子坐标取值范围
iniRange = r0Prm* (0.5+ (3*pointsNum/ (4* pi* sqrt(2)))^ (1/3) );

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
    while norm(randomPoint) > iniRange || minDis(1) < 0.9* r0Prm
        randomPoint = (2*rand(dimension,1)-1)*iniRange;
        minDis = distence(pointsList,randomPoint,dimension);
        %fprintf('太近了，重开\n')
    end
    pointsList((i-1)*dimension+1:i*dimension) = randomPoint;
end

fprintf('initial\n')

load('check_pointsList_mid.mat');

neighbor = neighborMatrix(pointsList,rcPrm,dimension,pointsNum);
eleDensity = electrondensity(pointsList,r0Prm,betaPrm,neighbor,dimension,pointsNum);
energy = potential(pointsList,r0Prm,E0Prm,Phi0Prm,alphaPrm,betaPrm,gammaPrm,neighbor,eleDensity,dimension,pointsNum);
gr = gradient(pointsList,r0Prm,E0Prm,Phi0Prm,alphaPrm,betaPrm,gammaPrm,neighbor,eleDensity,dimension,pointsNum);
dr = - gr;

% load('check_dr.mat')

%取样点
sample = 100;

%用数组记录grad,direction,energy,h,方向变化间隔
k = 1;
t = 1;
glist = gr;
dlist = dr;
elist = energy;
%搜索步长
% h = 1.2* max(abs(gr))/(4*sample);
h = 0.05617246/sample;
%fprintf('初始步长: %.4f\n',h)
% hlist = h;
% section = 1;
nantest = 0; %检测NaN报错
maxMove = 0.3;

drgrList = dr'* gr;
plist_hdr = h*dr;


for i = 0:sample
    pointsList = pointsList + plist_hdr;
    neighbor = neighborMatrix(pointsList,rcPrm,dimension,pointsNum);
    eleDensity = electrondensity(pointsList,r0Prm,betaPrm,neighbor,dimension,pointsNum);
    energy = potential(pointsList,r0Prm,E0Prm,Phi0Prm,alphaPrm,betaPrm,gammaPrm,neighbor,eleDensity,dimension,pointsNum);
    gr = gradient(pointsList,r0Prm,E0Prm,Phi0Prm,alphaPrm,betaPrm,gammaPrm,neighbor,eleDensity,dimension,pointsNum);

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





%所有近邻
function out = neighborMatrix(pointsList,rc,dim,num)
    out = zeros(num);
    for i = 1:num
        ri = pointsList((i-1)*dim+1:i*dim);
        for j = i+1:num
            if norm(ri-pointsList((j-1)*dim+1:j*dim)) <= rc
                out(i,j) = 1;
                out(j,i) = 1;
            end
        end
    end
end


%背景电子密度
function out = electrondensity(pointsList,r0,beta,neighbor,dim,num)
    out = zeros(num,1);
    for i = 1:num
        mid = 0;
        ri = pointsList((i-1)*dim+1:i*dim);
        rineighbor = neighbor(i,:);
        rineighbor = find(rineighbor);
        if ~isempty(rineighbor)
            for j = rineighbor
                mid = mid + exp( - beta* (norm(ri-pointsList((j-1)*dim+1:j*dim))/ r0 - 1) );
            end
            out(i) = mid/12;
        end
    end
end



%势能
function out = potential(pointsList,r0,E0,Phi0,alpha,beta,gamma,neighbor,eleDensity,dim,num)
    out = 0;

    %嵌入能
    for i = 1:num
        rhoi = eleDensity(i);
        if rhoi ~= 0
            out = out - E0* (1 - (alpha/beta)* log(rhoi))* rhoi^(alpha/beta) - 6* Phi0 * rhoi^(gamma/beta);
        end
    end
    
    %排斥能
    for i = 1:num-1
        ri = pointsList((i-1)*dim+1:i*dim);
        indices = find(neighbor(i,i+1:num)) + i;
        if ~isempty(indices)
            for j = indices 
                out = out + Phi0* exp( - gamma* (norm(ri-pointsList((j-1)*dim+1:j*dim))/ r0 - 1) );
            end
        end
    end
end

%梯度
function out = gradient(pointsList,r0,E0,Phi0,alpha,beta,gamma,neighbor,eleDensity,dim,num)
    out = zeros(dim*num,1);
    mid0 = zeros(dim,1);
        
    for i = 1:num

        rhoi =  eleDensity(i);
        ri = pointsList((i-1)*dim+1:i*dim);
        mid = mid0;

        if rhoi ~= 0
            rineighbor = neighbor(i,:);
            rineighbor = find(rineighbor);

            for j = rineighbor
                rhoj = eleDensity(j);
                rirj = ri - pointsList((j-1)*dim+1:j*dim);
                %嵌入能的导数
                mid = mid + (beta/ (12* r0))*...
                            ( - (alpha/beta)^2* E0* log(rhoi)* rhoi^ (alpha/beta- 1)+ 6* (gamma/beta)* Phi0* rhoi^ (gamma/beta- 1)...
                              - (alpha/beta)^2* E0* log(rhoj)* rhoj^ (alpha/beta- 1)+ 6* (gamma/beta)* Phi0* rhoj^ (gamma/beta- 1)...
                            )* exp( -beta* (norm(rirj)/r0 - 1) )* rirj/ norm(rirj) ;
                %排斥能的导数
                mid = mid - (gamma/r0)* Phi0* exp( -gamma* (norm(rirj)/r0 - 1) )* rirj/ norm(rirj);
            end
        end

        out((i-1)*dim+1:i*dim) = mid;
    end
end


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


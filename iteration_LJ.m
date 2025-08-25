%LJ势能

clear
set(groot , 'DefaultAxesFontName', '宋体', 'DefaultAxesFontSize', 15);  % 设置坐标轴标签的默认字体和大小
set(groot, 'defaultTextFontName', 'Times New Roman', 'defaultTextFontSize', 14);
% %画图参数
% set(0, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontSize', 20);  % 设置坐标轴标签的默认字体和大小
% set(0, 'DefaultTextFontSize', 20);  % 设置文本的默认字体大小
% set(0, 'DefaultLineLineWidth', 2);  % 设置默认线条的粗细
% set(0, 'DefaultLineMarkerSize', 12);  % 设置默认点的半径
% %set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',20,'DefaultTextFontSize',20,'DefaultLineLineWidth',2,'DefaultLineMarkerSize',15);

%系统编号
systemNum = 1;

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
        %fprintf('太近了，重开\n')
    end
    pointsList((i-1)*dimension+1:i*dimension) = randomPoint;
end

load('pointsList_1.mat');
save(['pointsList_',num2str(systemNum)],'pointsList')


energy = potential(epsilon,sigma,pointsList,dimension,pointsNum);
gr = gradient(epsilon,sigma,pointsList,dimension,pointsNum);%初始点梯度
dr = - gr;%初始搜索方向为负梯度方向

%用数组记录grad,direction,energy,h,方向变化间隔
k = 1;
t = 1;
glist = gr;
dlist = dr;
elist = energy;
%搜索步长
h = max(abs(gr))/100;
%fprintf('初始步长: %.4f\n',h)
%hlist = h;
%section = 1;
nantest = 0; %检测NaN报错
maxMove = 0.1;

%三阶插值近似求解的步长和精度
a = 0;
b = 0;
change = 0;%函数值交换用的中间值
accu = 0.1;
s = 0;
z = 0;
w = 0;

%循环次数
times = 0;
timeTotal = 0;

while norm(gr) > 0.005 && times < 1000
    times = times + 1;
    %fprintf('\n第 %d 个点,\n',times);

    %三阶插值
    t1 = 0;
    f1 = potential(epsilon,sigma,pointsList,dimension,pointsNum);
    g1 = gr'*dr;
    %fprintf('t1: %.8f, f1: %.4f, g1: %.4f, \n',t1,f1,g1);

    limit = 0;
    while (abs(g1) > accu && limit < 10)|| limit == 0
        limit = limit + 1;

        if g1 < 0
            h = abs(h);
        else
            h = - abs(h);
        end

        t2 = t1 + h;
        f2 = potential(epsilon,sigma,pointsList+t2*dr,dimension,pointsNum);
        g2 = gradient(epsilon,sigma,pointsList+t2*dr,dimension,pointsNum)'*dr;
        %fprintf('t2: %.8f, f2: %.4f, g2: %.4f,\n',t2,f2,g2)
        
        if abs(g2) < accu %&& abs(g2) < abs(g1)
            t1 = t2*2^sign(g1*g2);
            f1 = f2;
            g1 = g2;
            %fprintf('结束，点已满足\n')
            %fprintf('点： %.8f , f: %.6f , g: %.6f \n',t1,f1,g1);
            timeTotal = timeTotal + 1;
            elist = [elist,f1];
            %h = t1;
            %hlist = [hlist,h];
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

            fNew = potential(epsilon,sigma,pointsList+tNew*dr,dimension,pointsNum);
            gNew = gradient(epsilon,sigma,pointsList+tNew*dr,dimension,pointsNum)'*dr;
            %fprintf('求解出的新点：a1: %.8f, f1: %.4f, g1: %.4f, \n',tNew,fNew,gNew);

            if ( t1<tNew && tNew<t2 )||( t1>tNew && tNew>t2 )
                if gNew * g1 > 0 || abs(gNew) < abs(g1)/2
                    t1 = tNew;
                    f1 = fNew;
                    g1 = gNew;
                    %fprintf('取了新的点！\n')
                else
                    %fprintf('没有取新的点…\n')
                    h = h*1.5;%少收缩一些
                end
                
                if abs(gNew) > accu || limit == 1
                    h = h/6;
                    %fprintf('步长 h 太长了！！ %.5f\n',h);
                end
            else
                h = h/4;
                t1 = (t1+t2)/2;
                f1 = potential(epsilon,sigma,pointsList+t1*dr,dimension,pointsNum);
                g1 =  gradient(epsilon,sigma,pointsList+t1*dr,dimension,pointsNum)'*dr;
                %fprintf('三次插值出错，取中点\n')
            end
        elseif abs(g2) < 10*abs(g1)
            h = 2*h;
            t1 = t2;
            f1 = f2;
            g1 = g2;
            %fprintf('前进,而且步长 h 太短了…… %.5f\n',h);
        else
            h = h/6;
            %fprintf('梯度异常，返回\n');
        end
        timeTotal = timeTotal + 1;
        elist = [elist,f1];
        %hlist = [hlist,h];
    end

    if t1 == 0
        t1 = t1 + h;
    end

    %fprintf('总循环数:%d\n',timeTotal)

    h = abs(t1);
    k = k + 1;
    %section = [section,timeTotal];

    %限制最大位移
    %fprintf('不限制下最大的移动距离: %.4f\n',max(abs(t1*dr)));
    move = t1*dr;
    indicesMax = abs(move) > maxMove;
    move(indicesMax) = sign(move(indicesMax)) * maxMove;
    pointsList = pointsList + move;

%     %拉回远离的原子
%     avgDisList = averageDistence2List(pointsList,dimension,pointsNum);
% 
%     meanDis = mean(avgDisList);
%     stdDis = std(avgDisList);
% 
%     indicesFar = abs(avgDisList - meanDis) > 3 * stdDis;
%     farList = find(indicesFar);
% 
%     for i = 1:length(farList)
%         indixi = farList(i);
%         farPi = pointsList((indixi-1)*dimension+1:indixi*dimension);
%         targetPi = distence(pointsList,farPi,dimension);
%         targetPi = targetPi(2);
%         targetPi = pointsList((targetPi-1)*dimension+1:targetPi*dimension);
%         pointsList((indixi-1)*dimension+1:indixi*dimension) = targetPi + 1.15*sigma*(farPi-targetPi)/norm(farPi-targetPi);
%         %fprintf('移动了 %d 个原子\n',i)
%     end

    gr = gradient(epsilon,sigma,pointsList,dimension,pointsNum);
    if any(isnan(gr))
        %disp('gr = NaN, 退出\n');
        break;
    end
    glist = [glist,gr];

    g_k = gr;
    y_k1 = g_k - glist(:,k-1);
    d_k1 = dlist(:,k-1);

    beta_MHS = g_k'* y_k1/ (y_k1'* d_k1)*...
               (1- (g_k'* d_k1)^2/ (norm(g_k)^2* norm(d_k1)^2 ));
    beta_DY = norm(g_k)^2/ (y_k1'* d_k1);
    beta_N = max(0,min(beta_DY,beta_MHS));

    dr = -(1+ beta_N* (g_k'* d_k1)/ norm(g_k)^2)* g_k+ beta_N*d_k1;
    dlist = [dlist,dr];

end

energy = potential(epsilon,sigma,pointsList,dimension,pointsNum);
fprintf('最终点能量：%.4f,梯度模长：%.6f\n',energy,norm(gr))

recordList = [energy;times;timeTotal;nantest];
filename = sprintf('_%.4f_%d_%d_%d.mat', energy, times, timeTotal,nantest);

fileID = fopen('Result.txt','w');
fprintf(fileID,'%d_%.4f_%d_%d_%d\n', systemNum, energy, times, timeTotal,nantest);

save(['Result_',num2str(systemNum)],'pointsList','elist','nantest','energy','times','timeTotal')


% % 原始能量下降曲线
% figure('Position', [40, 680, 700, 600]);
% xlist= 1:length(elist);
% plot(xlist,elist,xlist,elist,'r.',section,elist(section),'kx','MarkerSize',11)
% xlabel('迭代次数'), ylabel('体系能量');
% title('原始能量下降曲线')
% ylim([-550,-300])

% 能量下降曲线
figure('Position', [40, 680, 700, 600]);xlist= 1:length(elist);
elist2 = log(elist-min(elist)+0.001)/log(10);
plot(xlist,elist2,xlist,elist2,'r.')
xlabel(['迭代次数, 总点数：',num2str(times)]), ylabel('体系能量');
title(['能量下降曲线，总循环数：',num2str(timeTotal)])

% 画出所有点的坐标
figure('Position', [740, 680, 700, 600]);
midList = 1:3:n;
% plot3(pointsList(midList),pointsList(midList+1),pointsList(midList+2),'.');
scatter3(pointsList(midList), pointsList(midList+1), pointsList(midList+2), 20, 'b', 'filled'); 
title('所有点的位置图')

fprintf('finish!\n')

% % 步长的变化
% figure('Position', [50, 0, 2800, 600]);
% xlist= 1:length(hlist);
% normglist = sqrt(sum(glist.^2, 1))/pointsNum;
% plot(xlist,hlist,xlist,hlist,'r.',section,hlist(section),'kx',section,normglist,'bo','MarkerSize',11)
% xlabel('迭代次数'), ylabel('搜索步长');
% title('步长 h 变化曲线')
% ylim([0,0.1])

% writeVestaFile('output.vesta',pointsList,dimension,pointsNum)


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



%测试将函数值和导数的求解写成 function 代码

%体系维数+粒子数目
dimension = 3;
pointsNum = 5;

%力场参数
epsilon = 1;
sigma = 3;

%初始点取值范围
iniRange = 6;

% 生成随机点的坐标
pointsList = (2*rand(dimension*pointsNum,1)-1)*iniRange;

energy = potential(epsilon,sigma,pointsList,dimension,pointsNum);
dr = gradient(epsilon,sigma,pointsList,dimension,pointsNum);
disp(dr)
disp(energy)

midlist = 1:dimension:dimension*pointsNum;
plot3(pointsList(midlist),pointsList(midlist+1),pointsList(midlist+2),'r.')




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

function out = gradient(epsi,sigma,plist,dim,num)
out = zeros(dim*num,1);
for i = 0:num-1
    ri = plist(i*dim+1:(i+1)*dim);
    for k = 1:dim
        mid = 0;
        for j = 0:num-1
            if j~= i
                normrj = norm(ri-plist(j*dim+1:(j+1)*dim));
                mid = mid + 4*epsi* (plist(i*dim+k)-plist(j*dim+k))* (-12*sigma^12/normrj^14 + 6*sigma^6/normrj^8);
            end
        end
        out(i*dim+k) = mid;
    end
end
end


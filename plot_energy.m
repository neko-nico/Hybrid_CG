%绘制MLFFs的能量曲线

elist = load('elistSave_Malon.mat').elist;

xlist= 1:length(elist);
elist2 = log(elist-min(elist)+0.001)/log(10);
plot(xlist,elist2,xlist,elist2,'r.')
xlabel('迭代次数,搜索方向改变100次')
title('能量下降曲线，总循环数：184')


% plot(xlist,elist)




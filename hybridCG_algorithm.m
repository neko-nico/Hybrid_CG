%论文中的算法

dim = 4;
iteration = 5;

glist = rand( dim , iteration);
glist = [glist,rand(dim,1)];
disp(glist)

dlist = rand( dim , iteration);
dlist = [dlist,rand(dim,1)];
disp(dlist)

k = 2;

beta_MHS = glist(:,k)'* (glist(:,k)- glist(:,k-1))/ ((glist(:,k)- glist(:,k-1))'* dlist(:,k-1))*...
           (1- (glist(:,k)'* dlist(:,k-1))^2/ (norm(glist(:,k))^2* norm(dlist(:,k-1))^2 ));
beta_DY = norm(glist(:,k))^2/ ((glist(:,k)- glist(:,k-1))'* dlist(:,k-1));
disp(beta_MHS)
disp(beta_DY)

beta_N = max(0,min(beta_DY,beta_MHS));
disp(beta_N)

d_k = -(1+ beta_N* (glist(:,k)'* dlist(:,k-1))/ norm(glist(:,k))^2)* glist(:,k )+ beta_N*dlist(:,k-1);
disp(d_k)





% for beta_HS = gr
%     fprintf('beta_HS = %d\n', beta_HS);
% end




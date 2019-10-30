solver = 'AMG_Adapt';
prec = 'NULL';
epsilon = 1e-2;
nsize = 0;
% n = 6:12;
% for i = 1: length(n)
%     nsize = 2^n(i);
%     test_matrices('fruitfly-biogrid.mat',solver,prec,nsize,epsilon)
% end
% hard mat
test_matrices('boneS10.mat',solver,prec,nsize,epsilon)
test_matrices('fl2010.mat',solver,prec,nsize,epsilon)
test_matrices('com-dblp.mat',solver,prec,nsize,epsilon)
test_matrices('youtube.mat',solver,prec,nsize,epsilon)
test_matrices('web-Google_Lg_comp.mat',solver,prec,nsize,epsilon)
test_matrices('web-BerkStan_Lg_comp.mat',solver,prec,nsize,epsilon)
test_matrices('amazon0601_Lg_comp.mat',solver,prec,nsize,epsilon)
test_matrices('web-NotreDame_Lg_comp.mat',solver,prec,nsize,epsilon)
test_matrices('web-Stanford_Lg_comp.mat',solver,prec,nsize,epsilon)
test_matrices('email-EuAll_Lg_comp.mat',solver,prec,nsize,epsilon)
test_matrices('loc-gowalla_edges_Lg_comp.mat',solver,prec,nsize,epsilon)

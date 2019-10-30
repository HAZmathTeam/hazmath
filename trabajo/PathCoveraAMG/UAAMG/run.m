solver = 'AMG_path';
prec = 'NULL';
epsilon
n = 6:11;
for i = 1: length(n)
    nsize = 2^n(i);
    test_matrices('fruitfly-biogrid.mat',solver,prec,nsize,epsilon)
end

%test_matrices('boneS10.mat',solver,prec)
%test_matrices('fl2010.mat',solver,prec)
%test_matrices('com-dblp.mat',solver,prec)
%test_matrices('youtube.mat',solver,prec)
%test_matrices('web-Google_Lg_comp.mat',solver,prec)
%test_matrices('web-BerkStan_Lg_comp.mat',solver,prec)
%test_matrices('amazon0601_Lg_comp.mat',solver,prec)
%test_matrices('web-NotreDame_Lg_comp.mat',solver,prec)
%test_matrices('web-Stanford_Lg_comp.mat',solver,prec)
%test_matrices('email-EuAll_Lg_comp.mat',solver,prec)
%test_matrices('loc-gowalla_edges_Lg_comp.mat',solver,prec)
%test_matrices('Ecoli-biogrid.mat',solver,prec)
%test_matrices('HIV-biogrid.mat',solver,prec)
%test_matrices('rat-biogrid.mat',solver,prec)
%test_matrices('fruitfly-biogrid.mat',solver,prec)
%test_matrices('human-biogrid.mat',solver,prec)
%test_matrices('mouse-biogrid.mat',solver,prec)
%test_matrices('worm-biogrid.mat',solver,prec)
%test_matrices('yeast-biogrid.mat',solver,prec)
%test_matrices('mouse.mat',solver,prec)
%test_matrices('human.mat',solver,prec)
%test_matrices('worm.mat',solver,prec)
%test_matrices('yeast.mat',solver,prec)


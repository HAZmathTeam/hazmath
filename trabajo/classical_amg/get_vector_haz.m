function [a]=get_vector_haz(chara)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% b=get_matrix_haz('b')
    %% loads a vector written in hazmath format.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% @Xiaozhe Hu, Department of Mathematics, Tufts University
    %% @Ludmil Zikatanov, Department of Mathematics, Penn State
    %--------------------------------------------------------------
    q0='''';
    disp(['Loading vector ',chara,'..'])
        a=eval(['load(',q0,chara,'.dat',q0,')']);
    nrow=a(1,1); 
    a=a(2:nrow+1,1:1);
    return
end


function [A]=get_matrix_haz(charA,sym)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Z=get_matrix_haz('Z',sym)
    %% loads matrix written in hazmath format.
    %% if general sym=0; if sym is not zero,
    %% then symmetric matrix is assumed and it is symmetrized
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% @Xiaozhe Hu, Department of Mathematics, Tufts University
    %% @Ludmil Zikatanov, Department of Mathematics, Penn State
    %--------------------------------------------------------------
    q0='''';
    q1 = '"';
    disp(['Loading matrix ',charA,'..'])
    %%    A=eval(['load(',q0,charA,'.dat',q0,')']);
    A=eval(['load(',q0,charA,q0,')']);
    nrow=A(1,1); ncol=A(1,2); nonz=A(1,3); A=A(2:nonz+1,1:3);
    A=sparse(A(:,1)+1,A(:,2)+1,A(:,3),nrow,ncol);
    if(sym), A=0.5*(A+A');  end    
    return
end


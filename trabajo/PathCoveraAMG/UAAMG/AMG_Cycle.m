function [x] = AMG_Cycle(amgData, b, x, level, amgParam)
% Multigrid cycle
%
% @ Xiaozhe Hu, Tufts University
    
% parameters
max_level = amgData(1).max_level;

n_presmooth = amgParam.n_presmooth;
n_postsmooth = amgParam.n_postsmooth;
cycle_type = amgParam.cycle_type;
ILU_level = amgParam.ILU_level;

% coarsest level
if (level == max_level) 
        
    x = ( amgData(level).A + 1.0e-10*speye(length(b), length(b)) )\b;
    %x = forward_gs(amgData(level).A, b, x, amgData(level).DL, 100);
    %x = backward_gs(amgData(level).A, b, x, amgData(level).DL, 100);

else
    
    % presmoothing
    if level <= ILU_level
        x = ILU_smoother(amgData(level).A, b, x, amgData(level).IL, amgData(level).IU, n_presmooth);
    else
        x = forward_gs(amgData(level).A, b, x, amgData(level).DL, n_presmooth);
        %x = jacobi(amgData(level).A, b, x, amgData(level).D, n_presmooth);
        %x = sym_gs(amgData(level).A, b, x, amgData(level).DL, amgData(level).DU, n_presmooth);
        %x = kaczmarz(amgData(level).A, b, x, n_presmooth, 0, 0);
    end
    
    % compute residual
    r = b - amgData(level).A*x;
    
    % restriction
    r_c = amgData(level).R*r;
                
    % coarse grid correction
    e_c = zeros(size(amgData(level+1).A,1),1);
    
    switch cycle_type
    
        case 'V',  
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
        
        case 'W',
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            
        case 'nV',  
            for i = 1:amgParam.coarse_it
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            end  
            
        case 'K',
            if (level == max_level-1)
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            else 
                e_0 = zeros(size(amgData(level+1).A,1),1);
                switch amgParam.coarse_krylov_type
                    case 'CG',                    
                        e_c = gcg(amgData(level+1).A, r_c, 1e-12, amgParam.coarse_it, ...
                            @(r)AMG_Cycle(amgData, r, e_0, level+1, amgParam), 0);
                    otherwise,
                        e_c = Prec_FGMRES(amgData(level+1).A, r_c, e_c, [], ...
                            @(r)AMG_Cycle(amgData, r, e_0, level+1, amgParam),...
                            amgParam.coarse_it, amgParam.coarse_it, 1e-12, 0);
                end
            end
            
        otherwise,
            display('Wrong cycle type!! run V-cycle!!')
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
    end
            
    % prolongation
    x = x + amgData(level).P*e_c;
        
    % postsmoothing
    if level <= ILU_level
        x = ILU_smoother(amgData(level).A, b, x, amgData(level).IL, amgData(level).IU, n_postsmooth);
    else
        x = backward_gs(amgData(level).A, b, x, amgData(level).DU, n_postsmooth);
        %x = jacobi(amgData(level).A, b, x, amgData(level).D, n_postsmooth);
        %x = sym_gs(amgData(level).A, b, x, amgData(level).DL, amgData(level).DU, n_presmooth);
        %x = kaczmarz(amgData(level).A, b, x, n_postsmooth, 0, 0);
    end
        
end

end

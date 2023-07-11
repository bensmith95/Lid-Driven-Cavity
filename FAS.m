%% FAS Alg.
% Inputs:
% Psi - Approximation for SF    W - Approximation for Vort.
% Sp - RHS of SF equation       Sw - RHS of Vort. eqn
% BC - Bound. Cond. function    h - size of grid
% Re - Reynolds Number          a - under-relaxation parameter
% tol - tolerance

% Outputs:
% Psi - Solution for SF at grid size h
% W - Solution for Vort. at grid size h

function [Psi,W] = FAS(Psi,W,Sp,Sw,BC,h,Re,a,tol)
    %% First Residual
    [Pres,Wres] = residual(Psi,W,Sp,Sw,h,Re); r = max([abs(Pres);abs(Wres)],[],'all');
    str = fprintf('residual = %f, %.2f%% Complete',[NaN,0]);
    
    %% Iterate until desired accuarcy
    while r>tol
        fprintf(repmat('\b',1,str))
        str = fprintf('residual = %f, %.2f%% Complete',[r,tol/r*100]);

        % V-Cycle
        [Psi,W] = Vcycle(Psi,W,Sp,Sw,BC,h,Re,a);

        % Residual
        [Pres,Wres] = residual(Psi,W,Sp,Sw,h,Re); 
        r = max([abs(Pres);abs(Wres)],[],'all');
    end
    fprintf(repmat('\b',1,str))
end
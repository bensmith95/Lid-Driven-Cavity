%% FMG-FAS Alg. to solve the non-linear system
% Inputs: 
% Psi - Initial Stream Function (SF) guess
% W - Initial Vorticity guess 
% BC - function that contains the Boundary Condition information
% H - size of fine grid
% Re - Reynolds Number

% Outputs:
% Psi - solution field for SF
% W - solution field for Vorticity
% res - final residual of solution

function [Psi,W,res] = FMG(Psi,W,BC,H,Re)
    %% Pre-Condition and Restrict to Coarse Grid
    tol = 1e-5; % desired accuracy/tolerance 
    a = 0.9; % under-relaxation parameter
    h = H; str1 = fprintf('Restricting to Coarse Grid...\n');
    while h<0.5
        % Relax
        Sp = zeros(size(Psi)); Sw = Sp;
        [Psi,W] = Relax(Psi,W,Sp,Sw,BC,h,Re,12,a);
        
        % Restrict
        Psi = restrict(Psi); W = restrict(W); h = 2*h; tol = 4*tol;
        N = size(Psi); Nx = N(2)-2; Ny = N(1)-2; hi = 1/h;
        
        % Update BCs
        [Ut] = BC(h,Re); 
        % vorticity 
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(2:Ny+1,Nx+1) = -2*Psi(2:Ny+1,Nx)*hi^2;
        W(2,2:Nx+1) = -2*Psi(3,2:Nx+1)*hi^2;
        W(Ny+1,2:Nx+1) = -2*Psi(Ny,2:Nx+1)*hi^2 - 2*Ut*hi;
        % Update extrapolated points
        Psi(2:Ny+1,1) = Psi(2:Ny+1,3);
        Psi(2:Ny+1,Nx+2) = Psi(2:Ny+1,Nx);
        Psi(1,2:Nx+1) = Psi(3,2:Nx+1);
        Psi(Ny+2,2:Nx+1) = 2*h*Ut + Psi(Ny,2:Nx+1);
        W(2:Ny+1,1) = 3*W(2:Ny+1,2) - 3*W(2:Ny+1,3) + W(2:Ny+1,4);
        W(2:Ny+1,Nx+2) = 3*W(2:Ny+1,Nx+1) - 3*W(2:Ny+1,Nx) + W(2:Ny+1,Nx-1);
        W(1,2:Nx+1) = 3*W(2,2:Nx+1) - 3*W(3,2:Nx+1) + W(4,2:Nx+1);
        W(Ny+2,2:Nx+1) = 3*W(Ny+1,2:Nx+1) - 3*W(Ny,2:Nx+1) + W(Ny-1,2:Nx+1);
    end
    %% Solve on Coarse Grid
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving on %iX%i grid\n',[Nx,Ny]);
    Sp = zeros(size(Psi)); Sw = Sp;
    [Pres,Wres] = residual(Psi,W,Sp,Sw,h,Re); r = max([abs(Pres);abs(Wres)],[],'all');
    str2 = fprintf('residual = %f, %.2f%% Complete',[r,tol/r*100]);
    while r>1e-8
        [Psi,W] = Relax(Psi,W,Sp,Sw,BC,h,Re,30,0.95);
        [Pres,Wres] = residual(Psi,W,Sp,Sw,h,Re); 
        r = max([abs(Pres);abs(Wres)],[],'all');
    end
    fprintf(repmat('\b',1,str2))
    %% FMG Alg.
    while h>H
        % Interpolate 
        Psi = interp(Psi); W = interp(W); h = h/2; tol = tol/4;
        N = size(Psi); Nx = N(2)-2; Ny = N(1)-2; hi = 1/h;
        fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving on %iX%i grid\n',[Nx,Ny]);
        
        % Update BCs
        [Ut] = BC(h,Re); 
        % vorticity 
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(2:Ny+1,Nx+1) = -2*Psi(2:Ny+1,Nx)*hi^2;
        W(2,2:Nx+1) = -2*Psi(3,2:Nx+1)*hi^2;
        W(Ny+1,2:Nx+1) = -2*Psi(Ny,2:Nx+1)*hi^2 - 2*Ut*hi;
        % Update extrapolated points
        Psi(2:Ny+1,1) = Psi(2:Ny+1,3);
        Psi(2:Ny+1,Nx+2) = Psi(2:Ny+1,Nx);
        Psi(1,2:Nx+1) = Psi(3,2:Nx+1);
        Psi(Ny+2,2:Nx+1) = 2*h*Ut + Psi(Ny,2:Nx+1);
        W(2:Ny+1,1) = 3*W(2:Ny+1,2) - 3*W(2:Ny+1,3) + W(2:Ny+1,4);
        W(2:Ny+1,Nx+2) = 3*W(2:Ny+1,Nx+1) - 3*W(2:Ny+1,Nx) + W(2:Ny+1,Nx-1);
        W(1,2:Nx+1) = 3*W(2,2:Nx+1) - 3*W(3,2:Nx+1) + W(4,2:Nx+1);
        W(Ny+2,2:Nx+1) = 3*W(Ny+1,2:Nx+1) - 3*W(Ny,2:Nx+1) + W(Ny-1,2:Nx+1);
        
        % Post-Relax
        Sp = zeros(size(Psi)); Sw = Sp;
        [Psi,W] = Relax(Psi,W,Sp,Sw,BC,h,Re,15,a);
        
        % FAS Alg.
        [Psi,W] = FAS(Psi,W,Sp,Sw,BC,h,Re,a,tol);
    end
    %% Determine Residual
    fprintf(repmat('\b',1,str1));
    [Pres,Wres] = residual(Psi,W,Sp,Sw,h,Re); res = max([abs(Pres);abs(Wres)],[],'all');
end
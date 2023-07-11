%% Recursive V-Cycles to improve approximation of solutions
function [Psih,Wh] = Vcycle(Psih,Wh,Sph,Swh,BC,h,Re,a)
    %% Condition current grid
    % Relax
    [Psih,Wh] = Relax(Psih,Wh,Sph,Swh,BC,h,Re,15,a);

    % Compute Residual Errors
    [Presh,Wresh] = residual(Psih,Wh,Sph,Swh,h,Re);
 
    %% Restrict to Coarse grid
    Psi2h = restrict(Psih); W2h = restrict(Wh); 
    Pres2h = restrict(Presh); Wres2h = restrict(Wresh); 
    N2h = size(Psi2h); N2x = N2h(2)-2; N2y = N2h(1)-2; hii = 1/(2*h); 
    
    % Update BCs
    [Ut] = BC(2*h,Re);
    % vorticity 
    W2h(2:N2y+1,2) = -2*Psi2h(2:N2y+1,3)*hii^2;
    W2h(2:N2y+1,N2x+1) = -2*Psi2h(2:N2y+1,N2x)*hii^2;
    W2h(2,2:N2x+1) = -2*Psi2h(3,2:N2x+1)*hii^2;
    W2h(N2y+1,2:N2x+1) = -2*Psi2h(N2y,2:N2x+1)*hii^2 - 2*Ut*hii;
    % Update extrapolated points
    Psi2h(2:N2y+1,1) = Psi2h(2:N2y+1,3);
    Psi2h(2:N2y+1,N2x+2) = Psi2h(2:N2y+1,N2x);
    Psi2h(1,2:N2x+1) = Psi2h(3,2:N2x+1);
    Psi2h(N2y+2,2:N2x+1) = 2*h*Ut + Psi2h(N2y,2:N2x+1);
    W2h(2:N2y+1,1) = 3*W2h(2:N2y+1,2) - 3*W2h(2:N2y+1,3) + W2h(2:N2y+1,4);
    W2h(2:N2y+1,N2x+2) = 3*W2h(2:N2y+1,N2x+1) - 3*W2h(2:N2y+1,N2x) + W2h(2:N2y+1,N2x-1);
    W2h(1,2:N2x+1) = 3*W2h(2,2:N2x+1) - 3*W2h(3,2:N2x+1) + W2h(4,2:N2x+1);
    W2h(N2y+2,2:N2x+1) = 3*W2h(N2y+1,2:N2x+1) - 3*W2h(N2y,2:N2x+1) + W2h(N2y-1,2:N2x+1);
    
    % form new RHS
    [Sp2h,Sw2h] = RHS(Psi2h,W2h,Pres2h,Wres2h,2*h,Re);

    % stop recursion at smallest grid size, otherwise continue recursion
    if h>0.2
        [Ap2h,Aw2h] = Relax(Psi2h,W2h,Sp2h,Sw2h,BC,2*h,Re,30,0.95);
    else 
        [Ap2h,Aw2h] = Vcycle(Psi2h,W2h,Sp2h,Sw2h,BC,2*h,Re,a);
    end
    %% Interpolate upto fine grid
    % determine error
    Ep2h = Ap2h - Psi2h;
    Ew2h = Aw2h - W2h;
    % interpolate error
    Eph = interp(Ep2h);
    Ewh = interp(Ew2h);
    % correct
    Psih = Psih + Eph;
    Wh = Wh + Ewh;
    
    N = size(Psih); Nx = N(2)-2; Ny = N(1)-2; hi = 1/h;
    
    % Update BCs
    [Ut] = BC(h,Re);
    % vorticity 
    Wh(2:Ny+1,2) = -2*Psih(2:Ny+1,3)*hi^2;
    Wh(2:Ny+1,Nx+1) = -2*Psih(2:Ny+1,Nx)*hi^2;
    Wh(2,2:Nx+1) = -2*Psih(3,2:Nx+1)*hi^2;
    Wh(Ny+1,2:Nx+1) = -2*Psih(Ny,2:Nx+1)*hi^2 - 2*Ut*hi;
    % Update extrapolated points
    Psih(2:Ny+1,1) = Psih(2:Ny+1,3);
    Psih(2:Ny+1,Nx+2) = Psih(2:Ny+1,Nx);
    Psih(1,2:Nx+1) = Psih(3,2:Nx+1);
    Psih(Ny+2,2:Nx+1) = 2*h*Ut + Psih(Ny,2:Nx+1);
    Wh(2:Ny+1,1) = 3*Wh(2:Ny+1,2) - 3*Wh(2:Ny+1,3) + Wh(2:Ny+1,4);
    Wh(2:Ny+1,Nx+2) = 3*Wh(2:Ny+1,Nx+1) - 3*Wh(2:Ny+1,Nx) + Wh(2:Ny+1,Nx-1);
    Wh(1,2:Nx+1) = 3*Wh(2,2:Nx+1) - 3*Wh(3,2:Nx+1) + Wh(4,2:Nx+1);
    Wh(Ny+2,2:Nx+1) = 3*Wh(Ny+1,2:Nx+1) - 3*Wh(Ny,2:Nx+1) + Wh(Ny-1,2:Nx+1);

    % Post-Smoothing
    [Psih,Wh] = Relax(Psih,Wh,Sph,Swh,BC,h,Re,10,a);
end

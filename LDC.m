%% Solves the 2D Navier-Stokes for the Lid Driven Cavity 
% Inputs:
% Re - Reynolds Number      H - size of fine grid

function LDC(Re,H)
    %% Initialise
    str1 = fprintf('Initialising...\n');
    
    % set coords
    x = 0:H:1; Nx = length(x); 
    y = 0:H:1; Ny = length(y);

    % assign variable storage
    Psi = zeros(Ny+2,Nx+2); Omega = zeros(Ny+2,Nx+2); 

    % Obtain and update BC's
    hi = 1/H; [Ut] = BC(H,Re);
    % vorticity 
    Omega(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
    Omega(2:Ny+1,Nx+1) = -2*Psi(2:Ny+1,Nx)*hi^2;
    Omega(2,2:Nx+1) = -2*Psi(3,2:Nx+1)*hi^2;
    Omega(Ny+1,2:Nx+1) = -2*Psi(Ny,2:Nx+1)*hi^2 - 2*Ut*hi;
    % Update extrapolated points
    Psi(2:Ny+1,1) = Psi(2:Ny+1,3);
    Psi(2:Ny+1,Nx+2) = Psi(2:Ny+1,Nx);
    Psi(1,2:Nx+1) = Psi(3,2:Nx+1);
    Psi(Ny+2,2:Nx+1) = 2*H*Ut + Psi(Ny,2:Nx+1);
    Omega(2:Ny+1,1) = 3*Omega(2:Ny+1,2) - 3*Omega(2:Ny+1,3) + Omega(2:Ny+1,4);
    Omega(2:Ny+1,Nx+2) = 3*Omega(2:Ny+1,Nx+1) - 3*Omega(2:Ny+1,Nx) + Omega(2:Ny+1,Nx-1);
    Omega(1,2:Nx+1) = 3*Omega(2,2:Nx+1) - 3*Omega(3,2:Nx+1) + Omega(4,2:Nx+1);
    Omega(Ny+2,2:Nx+1) = 3*Omega(Ny+1,2:Nx+1) - 3*Omega(Ny,2:Nx+1) + Omega(Ny-1,2:Nx+1);

    %% FMG-FAS Alg.
    fprintf(repmat('\b',1,str1)); str1 = fprintf('FMG-FAS Algorithm...\n');
    [Psi,Omega,~] = FMG(Psi,Omega,@BC,H,Re); 
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solution Converged.\n');

    %% POST PROCESS
    str2 = fprintf('Post-Processing...\n');
    % Determine U via finite differences of SF
    U = zeros(Ny,Nx); 
    U(2:Ny-1,:) = (Psi(4:Ny+1,2:Nx+1)-Psi(2:Ny-1,2:Nx+1))*0.5*hi;
    U(:,1) = 0; U(:,Nx) = 0; U(1,:) = 0; U(Ny,:) = Ut;  
    % Determine V via finite differences of SF
    V = zeros(Ny,Nx);
    V(:,2:Nx-1) = -(Psi(2:Ny+1,4:Nx+1)-Psi(2:Ny+1,2:Nx-1))*0.5*hi;
    V(:,1) = 0; V(:,Nx) = 0;  V(1,:) = 0; V(Ny,:) = 0;
    
    % Solve Poisson Equation for Pressure
    P = Pressure(U,V,H); 

    % save data to file
    Vars{1} = U; Vars{2} = V; Vars{3} = Psi(2:Ny+1,2:Nx+1); 
    Vars{4} = Omega(2:Ny+1,2:Nx+1); Vars{5} = P;
    filename = ['Flows/LDC_Re=',num2str(Re),'.mat'];
    if exist('Flows','dir')==0; mkdir Flows; end
    save(filename, 'Vars', 'x', 'y')
    fprintf(repmat('\b',1,str2)); str2 = fprintf('Flow saved in %s\n', filename); pause(1)

    fprintf(repmat('\b',1,str2)); fprintf(repmat('\b',1,str1));
end
%% BOUNDARY CONDITIONS
function [Ut] = BC(H,Re)
   Ut = 1; 
end
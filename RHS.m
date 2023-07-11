%% Determines new Right Hand Side of equations for V-Cycles
% Inputs:
% Psi - Approximation for SF    W - Approximation for Vort.
% Pres - residual of SF eqn     Wres - residual of Vort. eqn
% h - size of grid              Re - Reynolds Number

% Outputs:
% Sp - new RHS of SF equation   Sw - new RHS of Vort. eqn

function [Sph,Swh] = RHS(Psi,W,Pres,Wres,h,Re)
    %% Initalise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2;
    Sph = zeros(N); Swh = zeros(N); 
    
    %% Determine A(v2h)
    Rei = 1/Re; hi = 1/h;
    for j = 3:Ny
        for i = 3:Nx
            dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi;
            dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;
            % SF terms
            if dPsidy>=0
                dWdx = (3*W(j,i)-4*W(j,i-1)+W(j,i-2))*0.5*hi; 
            else
                dWdx = -(3*W(j,i)-4*W(j,i+1)+W(j,i+2))*0.5*hi; 
            end
            if -dPsidx>=0
                dWdy = (3*W(j,i)-4*W(j-1,i)+W(j-2,i))*0.5*hi; 
            else
                dWdy = -(3*W(j,i)-4*W(j+1,i)+W(j+2,i))*0.5*hi; 
            end
            LapW = (W(j-1,i)+W(j,i-1)-4*W(j,i)+W(j,i+1)+W(j+1,i))*hi^2;

            % Vort. terms
            LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;

            % Determine new RHS
            Sph(j,i) = (dPsidy*dWdx - dPsidx*dWdy - Rei*LapW) + Pres(j,i);
            Swh(j,i) = (W(j,i) + LapPsi) + Wres(j,i);
        end
    end
end

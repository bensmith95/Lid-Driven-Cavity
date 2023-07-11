%% Gauss-Seidel-Newton block iterative scheme to relax/solve SF-Vort. equations
% Inputs:
% Psi - Approximation for SF    W - Approximation for Vort.
% Sp - RHS of SF equation       Sw - RHS of Vort. eqn
% BC - Bound. Cond. function    h - size of grid
% Re - Reynolds Number          n - # of iterations
% a - under-relaxation parameter

% Outputs:
% Psi - Updated approx. for SF  W - Updated approx. for Vort. 

function [Psi,W] = Relax(Psi,W,Sp,Sw,BC,h,Re,n,a)
    %% Initalise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2; 
    Rei = 1/Re; hi = 1/h; [Ut] = BC(h,Re);

    %% Iterate
    % direction of iterations follows overall direction of flow
    for k = 1:n
        for j = Ny:-1:3 
            for i = 3:Nx
                dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi;
                dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;

                % Internal Points
                A = zeros(2); b = zeros(2,1); % Jacobian

                % SF
                if dPsidy>=0
                    dWdx = (3*W(j,i)-4*W(j,i-1)+W(j,i-2))*0.5*hi; X = 3*dPsidy*0.5*hi;
                else
                    dWdx = -(3*W(j,i)-4*W(j,i+1)+W(j,i+2))*0.5*hi; X = -3*dPsidy*0.5*hi;
                end
                if -dPsidx>=0
                    dWdy = (3*W(j,i)-4*W(j-1,i)+W(j-2,i))*0.5*hi; Y = -3*dPsidx*0.5*hi;
                else
                    dWdy = -(3*W(j,i)-4*W(j+1,i)+W(j+2,i))*0.5*hi; Y = 3*dPsidx*0.5*hi;
                end
                LapW = (W(j-1,i)+W(j,i-1)-4*W(j,i)+W(j,i+1)+W(j+1,i))*hi^2;

                A(1,2) = X+Y+4*Rei*hi^2;
                b(1) = -(dPsidy*dWdx - dPsidx*dWdy - Rei*LapW - Sp(j,i));

                % Vort.
                LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;
                A(2,1) = -4*hi^2; A(2,2) = 1;
                b(2) = -(W(j,i) + LapPsi - Sw(j,i));

                % solve linear system & correct
                s = A\b;
                Psi(j,i) = Psi(j,i) + a*s(1); 
                W(j,i) = W(j,i) + a*s(2);

                % Update BC's 
                if ismember(i,[3,4])
                    W(j,2) = -2*Psi(j,3)*hi^2; 
                    W(j,1) = 3*W(j,2) - 3*W(j,3) + W(j,4); 
                    Psi(j,1) = Psi(j,3); 
                end
                if ismember(i,[Nx-1,Nx])
                    W(j,Nx+1) = -2*Psi(j,Nx)*hi^2;
                    W(j,Nx+2) = 3*W(j,Nx+1) - 3*W(j,Nx) + W(j,Nx-1);
                    Psi(j,Nx+2) = Psi(j,Nx); 
                end
                if ismember(j,[3,4]) 
                    W(2,i) = -2*Psi(3,i)*hi^2;
                    W(1,i) = 3*W(2,i) - 3*W(3,i) + W(4,i); 
                    Psi(1,i) = Psi(3,i); 
                end
                if ismember(j,[Ny-1,Ny])
                    W(Ny+1,i) = -2*Psi(Ny,i)*hi^2 - 2*Ut*hi;
                    W(Ny+2,i) = 3*W(Ny+1,i) - 3*W(Ny,i) + W(Ny-1,i); 
                    Psi(Ny+2,i) = 2*h*Ut + Psi(Ny,i); 
                end
            end
        end
        % Update BC's
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
end
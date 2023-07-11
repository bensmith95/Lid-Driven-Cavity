%% Determine Pressure distribution
% via the Pressure Laplace equation
function P = Pressure(U,V,h)
    %% Initialise
    [Ny,Nx] = size(U); hi = 1/h;
    
    % map to staggered grid
    Us1 = (U(1:Ny-1,:)+U(2:Ny,:))*0.5; Vs1 = (V(:,1:Nx-1)+V(:,2:Nx))*0.5;
    Vs2 = (V(1:Ny-1,:)+V(2:Ny,:))*0.5; Us2 = (U(:,1:Nx-1)+U(:,2:Nx))*0.5;
    Nx = Nx-1; Ny = Ny-1;

    % determine f
    dUdx = (Us1(:,2:end)-Us1(:,1:end-1))*0.5*hi; dUdy = (Us2(2:end,:)-Us2(1:end-1,:))*0.5*hi; 
    dVdx = (Vs2(:,2:end)-Vs2(:,1:end-1))*0.5*hi; dVdy = (Vs1(2:end,:)-Vs1(1:end-1,:))*0.5*hi;
    F = 2*(dUdx.*dVdy - dVdx.*dUdy)*h^2; 

    %% Laplace operator 
    e = ones(Nx,1); I = speye(Nx); Z = sparse(Nx,Nx); 
    D = spdiags([e,-4*e,e],-1:1,Nx,Nx); D(1,1) = -3; D(Nx,Nx) = -3;
    Db = spdiags([e,-3*e,e],-1:1,Nx,Nx); Db(1,1) = -2; Db(Nx,Nx) = -2; 
    A = [Db, I, repmat(Z,[1,Ny-2])];
    for j = 2:Ny-1
        A = [A; repmat(Z,[1,j-2]), I, D, I, repmat(Z,[1,Ny-3-j+2])];
    end
    A = [A; repmat(Z,[1,Ny-2]), I, Db];

    % reorganise B into a vector
    f = zeros(Nx*Ny,1);
    for j=1:Ny
       f(1+(j-1)*Nx:j*Nx) = F(j,:); 
    end
    
    % set corner pressure
    f(1) = 0; A(1,1) = 1; A(1,2:end) = 0;
    
    %% Solve linear system
    Pv = A\f;
    
    % reorganise solution into matrix form
    P = zeros(Ny,Nx);
    for j=1:Ny
       P(j,:) = Pv(1+(j-1)*Nx:j*Nx); 
    end
    
    %% Extrapolate to standard grid
    [x,y] = ndgrid(h/2:h:1-h/2); P = griddedInterpolant(x,y,P','spline');
    [x,y] = ndgrid(0:h:1); P = P(x,y)';
end

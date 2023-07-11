
clear all; clc; close all

% Reynolds number
Re = 700; 

% load data
file = ['Flows/LDC_Re=',num2str(Re),'.mat']; load(file)
U = Vars{1}; V = Vars{2}; SF = Vars{3}; vort = Vars{4}; P = Vars{5};

% display info
set(0,'units','pixels'); disp = get(0,'ScreenSize');
set(0,'defaulttextinterpreter','latex')
Lx = disp(3); Ly = disp(4); 

%% Plot flow variables
figure(1); t1 = tiledlayout(2,3); 
TT1 = ['Flow Variables: $Re =$ ',num2str(Re)]; title(t1,TT1,'interpreter','latex');

% plot vector field
nexttile(t1); fn = pcolor(x,y,sqrt(U.^2+V.^2)); hold on; 
fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
ii = 16; jj = 16; 
quiver(x(1:ii:end),y(1:jj:end),U(1:jj:end,1:ii:end),V(1:jj:end,1:ii:end),'color','k');
c = colorbar('location','northoutside'); colormap('jet');
set(c.Label,{'String','Interpreter'},{'$||\mathbf{u}||$','latex'});
ylabel('$y$','interpreter','latex','rotation',0);

% plot U component
nexttile(t1); fn = pcolor(x,y,U); 
fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
c = colorbar('location','northoutside'); colormap('jet');
set(c.Label,{'String','Interpreter'},{'$u(x,y)$','latex'}); 

% plot V component
nexttile(t1); fn = pcolor(x,y,V); 
fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
c = colorbar('location','northoutside'); colormap('jet');
set(c.Label,{'String','Interpreter'},{'$v(x,y)$','latex'});  

% plot streamfunction
nexttile(t1); contourf(x,y,SF); 
c = colorbar('location','northoutside'); colormap('jet');
set(c.Label,{'String','Interpreter'},{'$\psi(x,y)$','latex'}); 
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex','rotation',0);

% plot vorticity
nexttile(t1); contourf(x,y,vort); 
c = colorbar('location','northoutside'); colormap('jet');
set(c.Label,{'String','Interpreter'},{'$\omega(x,y)$','latex'}); 
xlabel('$x$','interpreter','latex');

% plot pressure
nexttile(t1); fn = pcolor(x,y,P); 
fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
c = colorbar('location','northoutside'); colormap('jet');
set(c.Label,{'String','Interpreter'},{'$p(x,y)$','latex'});
xlabel('$x$','interpreter','latex');

% settings
t1.TileSpacing = 'tight'; t1.Padding = 'compact'; 
set(gcf, 'Position',  [0.225*Lx, 0.075*Ly, 0.55*Lx, 0.825*Ly])

%% Ghia et al (1982) comparison
if ismember(Re,[100,400,1000])
    figure(2); t2 = tiledlayout(1,3);  
    TT2 = ['Ghia \textit{et al} (1982) comparison: $Re =$ ',num2str(Re)]; 
    title(t2,TT2,'interpreter','latex');
    
    % load data
    if Re==100; c=2; elseif Re==400; c=3; else; c=4; end
    D1 = readmatrix('GhiaData1982.xlsx','Sheet','u geom mid','Range','A3:D19');
    D2 = readmatrix('GhiaData1982.xlsx','Sheet','v geom mid','Range','A3:D19');
    D3 = readmatrix('GhiaData1982.xlsx','Sheet','vort along U=1','Range','A3:D17');
    D4 = readmatrix('GhiaData1982.xlsx','Sheet','primary vortex','Range','B2:E4'); 

    % plot U(x=0.5,y)
    [~,xm] = min(abs(x-0.5));
    nexttile(t2); plot(y,U(:,xm)); grid on; hold on; plot(D1(:,1),D1(:,c),'kx'); 
    xlabel('$y$'); title('$u(0.5,y)$') 

    % plot V(x,y=0.5)
    [~,ym] = min(abs(y-0.5));
    nexttile(t2); plot(x,V(ym,:)); grid on; hold on; plot(D2(:,1),D2(:,c),'kx');
    xlabel('$x$'); title('$v(x,0.5)$')

    % plot vorticity along moving edge (U=1)
    nexttile(t2); plot(x(2:end-1),vort(end,2:end-1)); grid on; hold on; plot(D3(:,1),-D3(:,c),'kx');
    xlabel('$x$'); title('$\omega(x,1)$')

    % settings
    t2.TileSpacing = 'compact'; t2.Padding = 'tight'; 
    set(gcf, 'Position',  [0.075*Lx, 0.225*Ly, 0.825*Lx, 0.55*Ly])
    
    % position of primary vortex
    figure(3); [ypv,xpv] = find(SF==min(SF,[],'all')); 
    RowNames = {'Ghia et al (1982)','Smith (2023)'};
    ColNames = {[char(968),' min'],char(969),'x','y'};
    SFmin = [D4(c-1,1); SF(ypv,xpv)]; Vortmin = [-D4(c-1,2); vort(ypv,xpv)];
    xcoord = [D4(c-1,3); x(xpv)]; ycoord = [D4(c-1,4); y(ypv)];
    T = table(SFmin,Vortmin,xcoord,ycoord,'RowNames',RowNames);
    uitable('Data', T{:,:}, 'ColumnName', ColNames, 'RowName',T.Properties.RowNames, ...
            'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
    title = ['Position of Primary Vortex: Re = ',num2str(Re)];
    uicontrol('Style', 'text', 'Position', [20 300 200 15], 'String', title);
     
end


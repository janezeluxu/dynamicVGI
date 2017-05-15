p1 = 2;
p2 = 1;
p3 = 3;

% [X,Y] = meshgrid(-1:0.05:0,-1:0.05:0);
% Z = (X+1).^p1.*(Y+1).^p2;
% %mesh(X,Y,Z)
% hold on
% 
% [X,Y] = meshgrid(-1:0.05:0,-1:0.05:0);
% Z = -2*X.*(X+1).*(Y + 1);
% %mesh(X,Y,Z)
% hold on
% 
% % [X,Y] = meshgrid(0:0.05:1,0:0.05:1);
% % Z = (1-X).^p3.*(Y).^p2;
% % mesh(X,Y,Z)
% % hold on
% % 
% % [X,Y] = meshgrid(0:0.05:1,0:0.05:1);
% % Z = (1-X).^p3.*(1-Y).^p2;
% % mesh(X,Y,Z)
% % hold on

[t1,t2] = meshgrid(0:0.05:1,0:0.05:1);
lam1 = t1;
lam2 = t2.*(1-t1);
lam3 = (1-t1).*(1-t2);

%N1 = lam1- 0.75*(4*lam1.^3.*lam3)-0.5*(6*lam1.^2.*lam3.^2)-0.25*(4*lam1.*lam3.^3);
N1 = lam1 - lam1.*lam2-3*lam3.*lam1.^3-3*lam3.^2.*lam1.^2-lam3.^3.*lam1;
N2 = lam2-lam2.*lam3 - lam1.*lam2;
N3 = lam3-lam2.*lam3- lam1.^3.*lam3-3*lam1.^2.*lam3.^2-3*lam1.*lam3.^3;
N1_act = lam1.^4;
N2_act = lam2.^2;
N3_act = lam3.^2;

x = [0,-1,0];
y = [1,0,0];
xc = lam1*x(1)+lam2*x(2)+lam3*x(3);
yc = lam1*y(1)+lam2*y(2)+lam3*y(3);

ze = 2*(lam1.^0).*(lam2.^1).*(lam3.^1);

figure (1)
mesh(xc,yc,N1);
hold on
mesh(xc,yc,N2);
hold on
mesh(xc,yc,N3);
%mesh(xc,yc,N1_act);
%map = [0., 0., 0.];
%colormap(map)

% in t space
z1 = t1.^p2.*(1-t2)+t1.^p3.*t2;
z2 =(1-t1).^p3.*t2.^p1;
z3 =(1-t1).^p2.*(1-t2).^p1;

%lamda space
zl1 = (lam1.^p3.*lam2+lam1.^p2.*lam3)./(1-lam1);
zl2 = (lam2.^p1.*lam3+lam2.^p3.*lam1)./(1-lam2);
zl3 = (lam3.^p2.*lam1+lam3.^p1.*lam2)./(1-lam3);

zz1 = (lam1).*((1-lam2).^(p3-1)).*((1-lam3).^(p2-1)); %vertex 1
zz2 = (lam2).*((1-lam3).^(p1-1)).*((1-lam1).^(p3-1)); %vertex 2
zz3 = (lam3).*((1-lam1).^(p2-1)).*((1-lam2).^(p1-1)); %vertex 3

figure (2)
 mesh(xc,yc,z1);
 hold on
 mesh(xc,yc,z2);
 hold on
 mesh(xc,yc,z3);
 hold on

% mesh(xc,yc,zl1);
% hold on
% mesh(xc,yc,zl2);
% hold on
% mesh(xc,yc,zl3);
% hold on
% 
% mesh(xc,yc,zz1);
% hold on
% mesh(xc,yc,zz2);
% hold on
% mesh(xc,yc,zz3);
% hold on

%map = [0., 0., 0.];
%colormap(map)
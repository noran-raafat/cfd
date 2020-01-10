clc
clear all 
close all
tic
%% airfoil coordinates NACA 0012

%X_airfoil=[1	0.998973	0.995895	0.99078	0.983647	0.974528	0.963458	0.950484	0.935659	0.919044	0.900707	0.880723	0.859175	0.83615	0.811745	0.786058	0.759196	0.731269	0.702392	0.672683	0.642264	0.61126	0.5798	0.548012	0.516026	0.483974	0.451988	0.4202	0.38874	0.357736	0.327317	0.297608	0.268731	0.240804	0.213942	0.188255	0.16385	0.140825	0.119277	0.099293	0.080956	0.064341	0.049516	0.036542	0.025472	0.016353	0.00922	0.004105	0.001027	0	0.001027	0.004105	0.00922	0.016353	0.025472	0.036542	0.049516	0.064341	0.080956	0.099293	0.119277	0.140825	0.16385	0.188255	0.213942	0.240804	0.268731	0.297608	0.327317	0.357736	0.38874	0.4202	0.451988	0.483974	0.516026	0.548012	0.5798	0.61126	0.642264	0.672683	0.702392	0.731269	0.759196	0.786058	0.811745	0.83615	0.859175	0.880723	0.900707	0.919044	0.935659	0.950484	0.963458	0.974528	0.983647	0.99078	0.995895	0.998973	1];
%Y_airfoil=[0	0.000149	0.000595	0.001332	0.002352	0.003642	0.005188	0.006973	0.008978	0.01118	0.01356	0.016093	0.018755	0.021524	0.024374	0.027281	0.030221	0.033169	0.036099	0.038986	0.041803	0.044522	0.047116	0.049556	0.051812	0.053854	0.055652	0.057178	0.058401	0.059296	0.059838	0.060006	0.059782	0.059152	0.058107	0.056644	0.054763	0.052472	0.049782	0.046708	0.043271	0.039493	0.0354	0.031017	0.026371	0.021488	0.016391	0.0111	0.005632	0	-0.005632	-0.0111	-0.016391	-0.021488	-0.026371	-0.031017	-0.0354	-0.039493	-0.043271	-0.046708	-0.049782	-0.052472	-0.054763	-0.056644	-0.058107	-0.059152	-0.059782	-0.060006	-0.059838	-0.059296	-0.058401	-0.057178	-0.055652	-0.053854	-0.051812	-0.049556	-0.047116	-0.044522	-0.041803	-0.038986	-0.036099	-0.033169	-0.030221	-0.027281	-0.024374	-0.021524	-0.018755	-0.016093	-0.01356	-0.01118	-0.008978	-0.006973	-0.005188	-0.003642	-0.002352	-0.001332	-0.000595	-0.000149	0];
%X_airfoil=X_airfoil-0.5;
filename_density_L = 'C:\Users\dell\Desktop\AERO\CFD\project\book1'; 
X_airfoil=xlsread(filename_density_L,'S1:S99');
Y_airfoil=xlsread(filename_density_L,'T1:T99');
%% required input
imax=length(X_airfoil);
jmax=length(Y_airfoil);
%% Grid Construction "Algebric O Grid"
theta=linspace(0,360,imax);
Ro=5;
xb=Ro.*cosd(theta);
yb=Ro.*sind(theta);
x_mat=zeros(imax,jmax);
y_mat=x_mat;
x_mat(:,1)=X_airfoil(:);
x_mat(:,imax)=xb(:);
y_mat(:,1)=Y_airfoil(:);
y_mat(:,jmax)=yb(:);

for j=2:jmax-1
    for i=1:imax
      dx=(x_mat(i,imax)-x_mat(i,1))./(imax-1);
      dy=(y_mat(i,imax)-y_mat(i,1))./(jmax-1);
      x_mat(i,j)=x_mat(i,j-1)'+dx;
      y_mat(i,j)=y_mat(i,j-1)'+dy;
      
    end
    
end
figure(2)
plot(x_mat,y_mat,'b');
hold on
plot(x_mat',y_mat','b');
% for i=1:imax
%     plot(x_mat(i,:),y_mat(i,:));
% end
% for j=1:jmax
% plot(x_mat(:,j),y_mat(:,j));
% end
%[X1,Y1]=meshgrid(x_mat,y_mat);
axis equal
figexp(2,'algebraicgrid.eps','eps')
%% Elliptic Grid Generation
nmax=50000;
d_zeta=1/(imax-1);
d_eta=1/(jmax-1);
rms=zeros(1,nmax);
for n=1:nmax
    x_old=x_mat;
    y_old=y_mat;
    
    for j=2:jmax-1
        for i=1:imax
            if i==1
                x_zeta=(x_mat(i+1,j)-x_mat(imax-1,j))/(2*d_zeta);
                x_eta=(x_mat(i,j+1)-x_mat(i,j-1))/(2*d_eta);
                y_zeta=(y_mat(i+1,j)-y_mat(imax-1,j))/(2*d_zeta);
                y_eta=(y_mat(i,j+1)-y_mat(i,j-1))/(2*d_eta);
                alpha=x_eta.^2+y_eta.^2;
                beta= x_zeta.*x_eta+y_zeta*y_eta;
                gamma=x_zeta.^2+y_zeta.^2;
                c1=alpha/(d_zeta^2);
                c2=gamma/(d_eta^2);
                c3=beta/(2*d_zeta*d_eta);
                x_mat(i,j)=(0.5/(c1+c2))*(c1*(x_mat(i+1,j)+x_mat(imax-1,j))+c2*(x_mat(i,j+1)+x_mat(i,j-1))-c3*(x_mat(i+1,j+1)-x_mat(i+1,j-1)+x_mat(imax-1,j-1)-x_mat(imax-1,j+1)));
                y_mat(i,j)=(0.5/(c1+c2))*(c1*(y_mat(i+1,j)+y_mat(imax-1,j))+c2*(y_mat(i,j+1)+y_mat(i,j-1))-c3*(y_mat(i+1,j+1)-y_mat(i+1,j-1)+y_mat(imax-1,j-1)-y_mat(imax-1,j+1)));
            elseif i>1 && i<imax 
                x_zeta=(x_mat(i+1,j)-x_mat(i-1,j))/(2*d_zeta);
                x_eta=(x_mat(i,j+1)-x_mat(i,j-1))/(2*d_eta);
                y_zeta=(y_mat(i+1,j)-y_mat(i-1,j))/(2*d_zeta);
                y_eta=(y_mat(i,j+1)-y_mat(i,j-1))/(2*d_eta);
                alpha=x_eta.^2+y_eta.^2;
                beta= x_zeta.*x_eta+y_zeta*y_eta;
                gamma=x_zeta.^2+y_zeta.^2;
                c1=alpha/d_zeta^2;
                c2=gamma/d_eta^2;
                c3=beta/(2*d_zeta*d_eta);
                x_mat(i,j)=(0.5/(c1+c2))*(c1*(x_mat(i+1,j)+x_mat(i-1,j))+c2*(x_mat(i,j+1)+x_mat(i,j-1))-c3*(x_mat(i+1,j+1)-x_mat(i+1,j-1)+x_mat(i-1,j-1)-x_mat(i-1,j+1)));
                y_mat(i,j)=(0.5/(c1+c2))*(c1*(y_mat(i+1,j)+y_mat(i-1,j))+c2*(y_mat(i,j+1)+y_mat(i,j-1))-c3*(y_mat(i+1,j+1)-y_mat(i+1,j-1)+y_mat(i-1,j-1)-y_mat(i-1,j+1)));
            else
                x_mat(i,j)= x_mat(1,j);
                y_mat(i,j)= y_mat(1,j);
        end
    end
    end
    rms(n)=sqrt(((sum(sum(x_mat-x_old)))^2+(sum(sum(y_mat-y_old)))^2)/(imax*jmax));
end
figure(3)
plot(x_mat,y_mat,'color','black');
hold on
plot(x_mat',y_mat','color','black');
axis equal
figexp(3,'ellipticgrid.eps','eps')
figure(4)
plot(log10(1:nmax),log10(rms))
figexp(4,'rmsgrid.eps','eps')
%% Solution Scheme 
% initiation 
epsi=zeros(imax,jmax);
c1=zeros(imax,jmax); c2=c1; c3=c1;
jacobian=zeros(imax,jmax);
epsi_eta=zeros(imax,jmax);
epsi_zeta=zeros(imax,jmax);

V_inf=87.644;
a=0;      
       %angle of attack deg
%BC
epsi(1,jmax)=0;
% Kutta condition
epsi(1:imax,1)=epsi(1,2) ;
rms=zeros(1,nmax);
for i=1:imax-1
    epsi(i+1,jmax)=epsi(i,jmax)+V_inf*cosd(a)*(y_mat(i+1,jmax)-y_mat(i,jmax))-V_inf*sind(a)*(x_mat(i+1,jmax)-x_mat(i,jmax));
end
for n=1:nmax
    epsi_old=epsi;
    for j=2:jmax-1
        for i=1:imax
            if i==1
                x_zeta(i,j)=(x_mat(i+1,j)-x_mat(imax-1,j))./(2*d_zeta);
                x_eta(i,j)=(x_mat(i,j+1)-x_mat(i,j-1))./(2*d_eta);
                y_zeta(i,j)=(y_mat(i+1,j)-y_mat(imax-1,j))./(2*d_zeta);
                y_eta(i,j)=(y_mat(i,j+1)-y_mat(i,j-1))./(2*d_eta);
                alpha(i,j)=x_eta(i,j).^2+y_eta(i,j).^2;
                beta(i,j)= x_zeta(i,j).*x_eta(i,j)+y_zeta(i,j).*y_eta(i,j);
                gamma(i,j)=x_zeta(i,j).^2+y_zeta(i,j).^2;
                jacobian(i,j)= x_zeta(i,j).*y_eta(i,j)-x_eta(i,j).*y_zeta(i,j);
                c1(i,j)=alpha(i,j)./(d_zeta^2);
                c2(i,j)=gamma(i,j)./(d_eta^2);
                c3(i,j)=beta(i,j)./(2*d_zeta.*d_eta);
                epsi(i,j)=(0.5/(c1(i,j)+c2(i,j))).*(c1(i,j).*(epsi(i+1,j)+epsi(imax-1,j))+c2(i,j).*(epsi(i,j+1)+epsi(i,j-1))-c3(i,j).*(epsi(i+1,j+1)-epsi(i+1,j-1)+epsi(imax-1,j-1)-epsi(imax-1,j+1)));
            elseif i>1 && i<imax 
                x_zeta(i,j)=(x_mat(i+1,j)-x_mat(i-1,j))./(2*d_zeta);
                x_eta(i,j)=(x_mat(i,j+1)-x_mat(i,j-1))./(2*d_eta);
                y_zeta(i,j)=(y_mat(i+1,j)-y_mat(i-1,j))./(2*d_zeta);
                y_eta(i,j)=(y_mat(i,j+1)-y_mat(i,j-1))./(2*d_eta);
                alpha(i,j)=x_eta(i,j).^2+y_eta(i,j).^2;
                beta(i,j)= x_zeta(i,j).*x_eta(i,j)+y_zeta(i,j).*y_eta(i,j);
                gamma(i,j)=x_zeta(i,j).^2+y_zeta(i,j).^2;
                jacobian(i,j)= x_zeta(i,j).*y_eta(i,j)-x_eta(i,j).*y_zeta(i,j);
                c1(i,j)=alpha(i,j)./(d_zeta^2);
                c2(i,j)=gamma(i,j)./(d_eta^2);
                c3(i,j)=beta(i,j)./(2*d_zeta.*d_eta);
                epsi(i,j)=(0.5/(c1(i,j)+c2(i,j))).*(c1(i,j).*(epsi(i+1,j)+epsi(i-1,j))+c2(i,j).*(epsi(i,j+1)+epsi(i,j-1))-c3(i,j).*(epsi(i+1,j+1)-epsi(i+1,j-1)+epsi(i-1,j-1)-epsi(i-1,j+1)));
            else 
        x_zeta(i,j)=x_zeta(1,j);
        x_eta(i,j)=x_eta(1,j);
        y_zeta(i,j)=y_zeta(1,j);
        y_eta(i,j)=y_eta(1,j);
        jacobian(i,j)=jacobian(1,j);
        epsi(i,j)=epsi(1,j);
            end
    end
end
epsi(:,1)=epsi(1,2);%Kutta condition
rms(n)=sqrt(((max(max(epsi-epsi_old)))^2/(imax*jmax)));
end
figure(9)
plot(log10(1:nmax),log10(rms))
figexp(9,'rmspsi.eps','eps')
% for j=1,jmax
     j=1;
        for i=1:imax
    if i==1
                x_zeta(i,j)=(x_mat(i+1,j)-x_mat(imax-1,j))./(2*d_zeta);
                x_eta(i,j)=(-x_mat(i,j+2)+4*x_mat(i,j+1)-3*x_mat(i,j))./(2*d_eta);
                y_zeta(i,j)=(y_mat(i+1,j)-y_mat(imax-1,j))./(2*d_zeta);
                y_eta(i,j)=(-y_mat(i,j+2)+4*y_mat(i,j+1)-3*y_mat(i,j))./(2*d_eta);
                jacobian(i,j)=(x_zeta(i,j).*y_eta(i,j)-x_eta(i,j).*y_zeta(i,j));
               epsi_eta(i,j)=(-epsi(i,j+2)+4*epsi(i,j+1)-3*epsi(i,j))./(2*d_eta); 
               epsi_zeta(i,j)=(epsi(i+1,j)-epsi(imax-1,j))./(2*d_zeta);
            elseif i>1 && i<imax 
                x_zeta(i,j)=(x_mat(i+1,j)-x_mat(i-1,j))./(2*d_zeta);
                x_eta(i,j)=(-x_mat(i,j+2)+4*x_mat(i,j+1)-3*x_mat(i,j))./(2*d_eta);
                y_zeta(i,j)=(y_mat(i+1,j)-y_mat(i-1,j))./(2*d_zeta);
                y_eta(i,j)=(-y_mat(i,j+2)+4*y_mat(i,j+1)-3*y_mat(i,j))./(2*d_eta);
               jacobian(i,j)=(x_zeta(i,j).*y_eta(i,j)-x_eta(i,j).*y_zeta(i,j));
               epsi_eta(i,j)=(-epsi(i,j+2)+4*epsi(i,j+1)-3*epsi(i,j))./(2*d_eta); 
               epsi_zeta(i,j)=(epsi(i+1,j)-epsi(i-1,j))./(2*d_zeta);
    else
        x_zeta(i,j)=x_zeta(1,j);
        x_eta(i,j)=x_eta(1,j);
        y_zeta(i,j)=y_zeta(1,j);
        y_eta(i,j)=y_eta(1,j);
        jacobian(i,j)=jacobian(1,j);
        epsi_eta(i,j)=epsi_eta(1,j);
        epsi_zeta(i,j)=epsi_zeta(1,j);
    end
        end

    j=jmax;
        for i=1:imax
        
    if i==1
                x_zeta(i,j)=(x_mat(i+1,j)-x_mat(imax-1,j))./(2*d_zeta);
                x_eta(i,j)=(x_mat(i,j-2)-4*x_mat(i,j-1)+3*x_mat(i,j))./(2*d_eta);
                y_zeta(i,j)=(y_mat(i+1,j)-y_mat(imax-1,j))./(2*d_zeta);
                y_eta(i,j)=(y_mat(i,j-2)-4*y_mat(i,j-1)+3*y_mat(i,j))./(2*d_eta);
                jacobian(i,j)=(x_zeta(i,j).*y_eta(i,j)-x_eta(i,j).*y_zeta(i,j));
               epsi_eta(i,j)=(epsi(i,j-2)-4*epsi(i,j-1)+3*epsi(i,j))./(2*d_eta); 
               epsi_zeta(i,j)=(epsi(i+1,j)-epsi(imax-1,j))./(2*d_zeta); 
            elseif i>1 && i<imax 
                x_zeta(i,j)=(x_mat(i+1,j)-x_mat(i-1,j))./(2*d_zeta);
                x_eta(i,j)=(x_mat(i,j-2)-4*x_mat(i,j-1)+3*x_mat(i,j))./(2*d_eta);
                y_zeta(i,j)=(y_mat(i+1,j)-y_mat(i-1,j))./(2*d_zeta);
                y_eta(i,j)=(y_mat(i,j-2)-4*y_mat(i,j-1)+3*y_mat(i,j))./(2*d_eta);
               jacobian(i,j)=(x_zeta(i,j).*y_eta(i,j)-x_eta(i,j).*y_zeta(i,j));
               epsi_eta(i,j)=(epsi(i,j-2)-4*epsi(i,j-1)+3*epsi(i,j))./(2*d_eta); 
               epsi_zeta(i,j)=(epsi(i+1,j)-epsi(i-1,j))./(2*d_zeta);
    elseif i==imax
        x_zeta(i,j)=x_zeta(1,j);
        x_eta(i,j)=x_eta(1,j);
        y_zeta(i,j)=y_zeta(1,j);
        y_eta(i,j)=y_eta(1,j);
        jacobian(i,j)=jacobian(1,j);
        epsi_eta(i,j)=epsi_eta(1,j);
        epsi_zeta(i,j)=epsi_zeta(1,j);
    end
        end
for j=2:jmax-1
    for i=1:imax
        if i==1
        epsi_eta(i,j)=(epsi(i,j+1)-epsi(i,j-1))./(2*d_eta);
        epsi_zeta(i,j)=(epsi(i+1,j)-epsi(imax-1,j))./(2*d_zeta);
        elseif i>1 && i<imax
        epsi_eta(i,j)=(epsi(i,j+1)-epsi(i,j-1))./(2*d_eta);
        epsi_zeta(i,j)=(epsi(i+1,j)-epsi(i-1,j))./(2*d_zeta);
        else
        epsi_eta(i,j)=epsi_eta(1,j);
        epsi_zeta(i,j)=epsi_zeta(1,j);
        end
        
    end
end
figure(6)
hold on
plot(x_mat(:,1),y_mat(:,1))
contour(x_mat,y_mat,epsi)
axis equal
figexp(6,'streamlines.eps','eps')
u=zeros(imax,jmax);
v=zeros(imax,jmax);
cp=zeros(imax,jmax);
for j=1:jmax
    for i=1:imax
        u(i,j)=(1./(jacobian(i,j))).*(epsi_zeta(i,j).*-x_eta(i,j)+epsi_eta(i,j).*x_zeta(i,j));
        v(i,j)=(1./(jacobian(i,j))).*(epsi_zeta(i,j).*-y_eta(i,j)+epsi_eta(i,j).*y_zeta(i,j));
        cp(i,j)=1-((u(i,j).^2+v(i,j).^2)./(V_inf.^2));
        V(i,j)=sqrt(u(i,j).^2+v(i,j).^2);
    end
end
figure(7)
plot(X_airfoil,-cp(:,1))
figexp(7,'cp.eps','eps')
cps=cp(:,1);
Cx=0;Cy=0;
for i=1:imax-1
    ex=(X_airfoil(i+1)-X_airfoil(i))./sqrt((X_airfoil(i+1)-X_airfoil(i)).^2+(Y_airfoil(i+1)-Y_airfoil(i)).^2);
    ey=(Y_airfoil(i+1)-Y_airfoil(i))./sqrt((X_airfoil(i+1)-X_airfoil(i)).^2+(Y_airfoil(i+1)-Y_airfoil(i)).^2);
    ds=sqrt((X_airfoil(i+1)-X_airfoil(i)).^2+(Y_airfoil(i+1)-Y_airfoil(i)).^2);
    nx=ey;
    ny=ex;
    Cx=Cx-0.5*(cps(i+1)+cps(i)).*nx.*ds;
    Cy=Cy+0.5*(cps(i+1)+cps(i)).*ny.*ds; 
end
cl=Cy.*cosd(a)-Cx.*sind(a);
cd=Cy.*sind(a)+Cx.*cosd(a);
figure(8)
hold on
plot(X_airfoil,Y_airfoil)
contourf(x_mat,y_mat,V,'linecolor','none')
axis equal
figexp(1,'velocitycontour.eps','eps')
figure()
figure;
plot(X_airfoil,V(:,1)/V_inf,'LineWidth',1);
grid on;
toc
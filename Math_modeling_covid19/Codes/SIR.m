function [t,X]=SIR()
nu=1/3;
beta=1/2;
eps=(1.25)*(10^(-6));

S0=1-eps; 
I0=eps;   
R0=0;    
Ti=0;   
Tf=100;
dt=1;
N=(Tf-Ti)/dt;
X=zeros(3,1);     
dX=zeros(3,1);   
t_tab=zeros(N+1,1);   
X_tab=zeros(3,N+1);  
t_tab(1,1)=Ti; 
X_tab(1,1)=S0;  
X_tab(2,1)=I0;  
X_tab(3,1)=R0;  
P_tab=zeros(N+1,1);
P_tab(1,1)=nu*I0;

for i=1:N
    t = t_tab(i,1);
    S = X_tab(1,i);
    I = X_tab(2,i);
    R = X_tab(3,i);
    X(1,1)=S; 
    X(2,1)=I;
    X(3,1)=R;
    dS = -beta*S*I;
    dI = beta*S*I-nu*I;
    dR = nu*I;
    dX(1,1)=dS;
    dX(2,1)=dI;
    dX(3,1)=dR;
    Xnew = X+dt*dX;
    X_tab(:,i+1)=Xnew;
    t_tab(i+1,1)=t+dt;
    P_tab(i+1,1)=nu*X_tab(2,i+1);  
end

figure(1)
plot(t_tab(:,1),X_tab(1,:),'LineWidth',2)
hold on 
plot(t_tab(:,1),X_tab(2,:),'LineWidth',2)
hold on 
plot(t_tab(:,1),X_tab(3,:),'LineWidth',2)
hold on
legend({'Susceptibles','Infectés','Rétablis'},'Location','northeast')

Int1=cumtrapz(t_tab(:,1),P_tab(:,1));
Int2 =trapz(t_tab(:,1),P_tab(:,1));

figure(2)
plot(t_tab(:,1),Int1,'LineWidth',2)
figure(3)
plot(t_tab(:,1),Int2,'LineWidth',2)

end
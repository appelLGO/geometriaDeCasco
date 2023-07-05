clear all, close all, clc

Cx1       = linspace(0.5,0.6244,5)';
Cx2       = linspace(0.6244,0.99,10)';
N_points  = 5001;

W         = 1;
K         = 0;
T         = 1;

slim_z    = zeros(N_points,length(Cx1));
slim_y    = zeros(N_points,length(Cx1));

for n = 1:length(Cx1)
  slim_z(:,n) = slimSectionIt(W,K,T,Cx1(n),N_points)(:,1);
  slim_y(:,n) = slimSectionIt(W,K,T,Cx1(n),N_points)(:,2);
endfor  

bunt_z    = zeros(N_points,length(Cx2));
bunt_y    = zeros(N_points,length(Cx2));

for n = 1:length(Cx2)
  bunt_z(:,n) = buntSectionIt(W,K,T,Cx2(n),N_points)(:,1);
  bunt_y(:,n) = buntSectionIt(W,K,T,Cx2(n),N_points)(:,2);
endfor

figure
plot(slim_y,slim_z,'b','linewidth',1.5,...
     bunt_y,bunt_z,'r','linewidth',1.5,...
     slim_y,zeros(N_points,1),':k',...
     slim_y,zeros(N_points,1) +1,':k',...
     zeros(N_points,1),slim_z,':k',...
     zeros(N_points,1) +1,slim_z,':k')
axis([-0.05 1.05 -0.05 1.05],'equal')
xlabel('\phi')
ylabel('\eta')
grid on

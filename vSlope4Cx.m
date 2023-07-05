clear all, close all, clc

Cx        = [0.5:0.001:1]';
Cx_switch = 0.6244;

W         = 1;
K         = 0;
T         = 1;
N_points  = 1e2;

for n = 1:length(Cx)
  if Cx(n) > Cx_switch
    z = buntSectionIt(W,K,T,Cx(n),N_points)(:,1);
    y = buntSectionIt(W,K,T,Cx(n),N_points)(:,2);
  elseif Cx(n) <= Cx_switch
    z = slimSectionIt(W,K,T,Cx(n),N_points)(:,1);
    y = slimSectionIt(W,K,T,Cx(n),N_points)(:,2);
  endif

  y1 = y(length(y));
  y0 = y(length(y) -1);
  z1 = z(length(z));
  z0 = z(length(z) -1);
  
  v_slope(n) = (y1 -y0)/(z1 -z0);
endfor

v_angle = atan(v_slope);

figure
  subplot(2,1,1)
    plot(Cx, v_slope,...
         zeros(length(Cx)) + Cx_switch,1.05*v_slope,':k') 
    grid on
    text(Cx_switch +0.01,0.022,num2str(Cx_switch))
    axis([0.5 1 0 1.05*max(v_slope)])
    xlabel('C_X')
    ylabel('\phi''_1')

  subplot(2,1,2)
    plot(Cx, 180/pi*v_angle,...
         zeros(length(Cx)) + Cx_switch,1.1*180/pi*v_angle,':k') 
    grid on
    text(Cx_switch +0.01,1.1,num2str(Cx_switch))
    axis([0.5 1 0 1.1*max(180/pi*v_angle)])
    xlabel('C_X')
    ylabel('\theta_1 [deg]')
    

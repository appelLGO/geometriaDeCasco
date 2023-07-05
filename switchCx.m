clear all, close all, clc

Cx = [0.5:0.025:0.75]';
N_points = 50;

z = linspace(0,1,N_points);
diff2_mean = zeros(length(Cx),1);
diff2_max = zeros(length(Cx),1);

for n = 1:length(Cx)
  slim           = slimSectionIt(1,0,1,Cx(n),N_points)(:,2);
  bunt           = buntSectionIt(1,0,1,Cx(n),N_points)(:,2);
  diff2(:,n)     = (bunt -slim).^2; 
  diff2_mean(n)  = mean(diff2(:,n));
  diff2_max(n)   = max(diff2(:,n)); 
  z_diff2_max(n) = z(find(diff2(:,n) == diff2_max(n))); 
endfor

min_diff2_mean = min(diff2_mean)
Cx_mean        = Cx(find(diff2_mean == min(diff2_mean))) 
min_diff2_max  = min(diff2_max)
Cx_max         = Cx(find(diff2_max == min(diff2_max)))
  
figure 1
subplot(1,2,1)
  plot(diff2,z,'k',...
       diff2_max,z_diff2_max,':or',...
       diff2_mean,z_diff2_max,':ob')
  grid on
  axis([0 0.077 0 1])
  title('\delta^2(\eta), \delta^2_{mean}, \delta^2_{max}, various C_X')
  ylabel('\eta')
  xlabel('\delta^2 (\eta)')
  text(diff2_max(1),z_diff2_max(1) +0.02,'Cx = 0.5')
  text(diff2_max(length(diff2_max)),z_diff2_max(length(z_diff2_max)) -0.02,'Cx = 0.75') 
 
  
subplot(1,2,2)
  plot(Cx,diff2_mean,'-b',... 
       Cx,diff2_max,'-r',...  
       Cx_mean, min_diff2_mean,'ob',...
       Cx_max, min_diff2_max,'or')
  grid on
  axis([0.49 0.76 -0.009 0.07])
  title('\delta^2_{mean}, \delta^2_{max}')
  ylabel('\delta^2_{mean}, \delta^2_{max} (C_X)')
  xlabel('C_X')
  legend('\delta^2_{mean}','\delta^2_{max}')
  text(Cx_mean, min_diff2_mean -0.002,'Cx = 0.6522')
  text(Cx_max, min_diff2_max +0.002,'Cx = 0.6244')


 
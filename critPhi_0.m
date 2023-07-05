clear all, close all, clc

eta = [0:1e-3:1]';

E       = 1;
count   = 0;
dPdn_u  = 5;
dPdn_l  = 0; 

while (E > 1e-6)&&(count < 100)
  dPdn  = 0.5*(dPdn_u +dPdn_l);
  Phii  = spline([0,1],[dPdn,0,1,0],eta);
  
  if max(Phii)>1
    dPdn_u = dPdn;
  else
    dPdn_l = dPdn;
  endif  
  
  E      = dPdn_u -dPdn_l;
  count += 1;
  
endwhile  
  

dPdn_0 = [0:0.5:6]';

figure
for n = 1:length(dPdn_0)
  Phi(:,n) = spline([0,1],[dPdn_0(n),0,1,0],eta);
  if max(Phi(:,n)) <= 1
    plot (eta,Phi(:,n),'b')
  elseif max(Phi(:,n)) > 1
    plot (eta,Phi(:,n),'r')
  endif
  hold on
endfor

grid on
xlabel('\eta')
ylabel('\phi (\eta)')
title('\phi (\eta, \phi^{,}_{0,n})')
hold on

% Curva p/ derivada crítica (dPdn = 3)
  hold on
  Phi_limit = spline([0,1],[dPdn,0,1,0],eta);

  plot(eta,Phi_limit,'k','linewidth',2)
  axis('equal')
  axis([-0.5 1.5 -0.5 1.5])



% Contorno do Retângulo 1x1
  hold on
  zero_vec      = zeros(length(eta),1);
  zero2one_vec  = linspace(0,1,length(eta));
  plot(eta,zero_vec,':k','linewidth',1.5, ...
       eta,zero_vec +1,':k','linewidth',1.5, ...
       zero_vec,zero2one_vec,':k','linewidth',1.5, ...
       zero_vec +1,zero2one_vec,':k','linewidth',1.5)



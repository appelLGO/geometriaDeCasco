clear all, close all, clc

x = [pi-1:0.001:pi+1]';
y = sqrt(x).*tan(x);

% Método da Bisseção

maxIt = 8;
xl = zeros(maxIt,1);
xu = zeros(maxIt,1);
xm = zeros(maxIt,1);

xl(1) = 2.5;
xu(1) = 4;

for n = 1:maxIt
  xm(n) = mean([xu(n),xl(n)]);
  yl    = sqrt(xl(n)).*tan(xl(n));
  ym    = sqrt(xm(n)).*tan(xm(n));
  
  if yl*ym > 0
    xu(n+1) = xu(n);
    xl(n+1) = xm(n);
  elseif yl*ym < 0 ;
    xu(n+1) = xm(n);
    xl(n+1) = xl(n);
  endif

endfor

yh = linspace(2.5,-1.5,maxIt)';

figure 
 plot(x,y,'linewidth',1.5,...
      xu(1:end-1),yh,'--<r',...
      xl(1:end-1),yh,'-->b',...
      xm,yh,'--ok',...
      zeros(maxIt,1)+pi,1.0*yh,'k')
  grid on
  axis([min(x)-0.25 max(x)+0.25 min(y)-0.1 max(y)+0.1]) 
  xlabel('x','fontsize',13)
  ylabel('y','fontsize',13)
  
for n = 1:maxIt
    hold on
    plot(linspace(2.3,xm(n)-0.05,10)',zeros(10,1)+yh(n),':k')
    text(2.3,yh(n)+0.12,strcat('n =', num2str(n)))
endfor    

  
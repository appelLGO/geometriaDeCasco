clear all, close all, clc

N = 2.^[-3:0.5:4]'; 
theta = [0:0.005*2*pi:pi/2]';

figure
hold on
xlabel('\eta')
ylabel('\phi')

for n = 1:length(N)
  R(:,n) = (cos(theta).^N(n) +sin(theta).^N(n)).^(-1/N(n));
  x(:,n) = real(R(:,n).*e.^(i*theta));
  y(:,n) = imag(R(:,n).*e.^(i*theta));

  plot(x(:,n),y(:,n),':k',...
       x(:,n),-y(:,n),':k',...
       -x(:,n),y(:,n),':k',...
       -x(:,n),-y(:,n),':k')
  grid on
  hold on
  axis([-1.5 1.5 -1.5 1.5],'equal')
  
  if N(n) >= 2
    plot(x(:,n),-y(:,n),'k','linewidth',1.5)
  endif
endfor


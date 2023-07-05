clear all, close all, clc

n = 1 .+ 2.^[-4:0.5:5]';
s = linspace(0,3,length(n))';

N_points = 50;

for m = 1:length(s)
z_slim = slimSection(s(m),N_points)(:,1);
y_slim = slimSection(s(m),N_points)(:,2);

Cx_slim(m) = trapz(z_slim,y_slim);
endfor

for m = 1:length(n)
z_bunt = buntSection(n(m),N_points)(:,1);
y_bunt = buntSection(n(m),N_points)(:,2);

Cx_bunt(m) = trapz(z_bunt,y_bunt);
endfor

figure 1

subplot(1,2,1)
  plot(s,Cx_slim,'-ob',...
       s,zeros(length(s)) +0.5,':k',...
       s,zeros(length(s)) +0.75,':k')
  text([0.1,0.1],[0.47,0.78],{'C_X = 0.5','C_X = 0.75'}) 
  grid on
  axis([ -0.2 1.1*max(s) 0 1.09])
  xlabel('\phi''_0')
  ylabel('C_X')

subplot(1,2,2)
plot(n,Cx_bunt,'-or',...
     n,zeros(length(n)) +0.5,':k',...
     n,zeros(length(n)) +1,':k')
  text([0.1,0.1],[0.47,1.03],{'C_X = 0.5','C_X = 1'})
grid on
axis([ min(n)-3 1.1*max(n) 0 1.09])
xlabel('s')
ylabel('C_X')
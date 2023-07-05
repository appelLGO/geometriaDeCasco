clear all,close all, clc

% Rotina para Geração de Geometria de Casco por Curvas Paramétricas

% Grandezas Iniciais & Parâmetros para ajuste geométrico
    Lbp  = 12;      
    B    = 2.7;    
    T    = 1.3;    
    Cm   = 0.9;    

    kL_pmb      = 4/11;
    kX_pmb      = 0.45;
    
    k_SAC_aft   = 0.1;
    k_SAC_fwd   = 0.05;
    ks_SAC_aft  = 1.5;
    ks_SAC_fwd  = 0.05;
    
    k_K_aft     = 0.85;
    k_K_fwd     = 0.45; 
    ks_K_aft    = 2.75;   
    ks_K_fwd    = 2.5;
    
    k_W_aft     = 0.9;     
    k_W_fwd     = 0.1;     
    ks_W_aft    = 2.5;   
    ks_W_fwd    = 1.5; 

% ======================================================================================================    
% Corpo Médio Paralelo e Vetor de Coord. Longitudinal
    
    L_pmb = kL_pmb*Lbp;               % [m] Comprim. do Corpo Médio Paralelo  
    X_pmb  = kX_pmb*Lbp;              % [m] Coord. Long. do Centro do Corpo Médio Paralelo  
    
    X_aft = ceil(X_pmb - 0.5*L_pmb);  % [m] Coord. Long. do Extremo a Ré do Corpo Médio Paralelo
    X_fwd = floor(X_pmb + 0.5*L_pmb); % [m] Coord. Long. do Extremo a Vante do Corpo Médio Paralelo 
   
    x_aft  = [0:1:X_aft]';
    x_pmb0 = [X_aft:1:X_fwd]';
    x_pmb  = x_pmb0(2:(length(x_pmb0)-1));
    x_fwd  = [X_fwd:1:Lbp]';

  x = [x_aft; x_pmb; x_fwd]; 

% SAC
    SAC_aft   = k_SAC_aft*B*T*Cm;                               % [m2]    Área esp. de popa
    SAC_fwd   = k_SAC_fwd*B*T*Cm;                               % [m2]    Área esp. de proa
    s_SAC_aft = ks_SAC_aft*((B*T*Cm -SAC_aft)/X_aft);           % [m2/m]  Deriv. SAC esp. popa
    s_SAC_fwd = ks_SAC_fwd*((SAC_fwd -B*T*Cm)/(Lbp -X_fwd));    % [m2/m]  Deriv. SAC esp. proa      
    
    % Corpo de Ré
    SAC_aft = spline([0, X_aft],[s_SAC_aft, SAC_aft, B*T*Cm, 0], x_aft);
    
    % Corpo de Vante
    SAC_fwd = spline([X_fwd, Lbp],[0, B*T*Cm, SAC_fwd, s_SAC_fwd], x_fwd);
   
    % Corpo Médio Paralelo 
    SAC_pmb = B*T*Cm .+ zeros(size(x_pmb));

  SAC = [SAC_aft;SAC_pmb;SAC_fwd];

% Curva da quilha
    K_aft       = k_K_aft*T;                      % [m]    Elevação da quilha À Ré 
    K_fwd       = k_K_fwd*T;                      % [m]    Elevação a Vante 
    s_aft_keel  = ks_K_aft*(-K_aft/X_aft);        % [m/m]  Deriv. Quilha esp. popa
    s_fwd_keel  = ks_K_fwd*(K_fwd/(Lbp -X_fwd));  % [m/m]  Deriv. Quilha esp. proa      
    
    % Corpo de Ré
    K_aft = spline([0 , X_aft],[s_aft_keel, K_aft, 0, 0], x_aft);
    
    % Corpo de Vante
    K_fwd = spline([X_fwd, Lbp],[0, 0, K_fwd, s_fwd_keel], x_fwd);
    
    % Corpo Médio Paralelo
    K_pmb = zeros(size(x_pmb));

  K = [K_aft;K_pmb;K_fwd];    
 
% Perfil de Plano de Linha D'água 
    W_aft   = k_W_aft*B/2;                         % [m]    Meia-boca a Ré (Propulsor)
    W_fwd   = k_W_fwd*B/2;                         % [m]    Meia-boca a Vante (Bulbo e Roda de Proa)
    s_W_aft = ks_W_aft*(B/2 -W_aft)/X_aft;         % [m/m]  Deriv. Quilha esp. popa
    s_W_fwd = ks_W_fwd*(W_fwd -B/2)/(Lbp -X_fwd);  % [m/m]  Deriv. Quilha esp. proa      
    
    % Corpo de Ré
    W_aft = spline([0, X_aft],[s_W_aft, W_aft, B/2, 0], x_aft);
    
    % Corpo de Vante
    W_fwd = spline([X_fwd, Lbp],[0, B/2, W_fwd, s_W_fwd], x_fwd);
   
    % Corpo Médio Paralelo  
    W_pmb = B/2 .+zeros(size(x_pmb));

  W = [W_aft;W_pmb;W_fwd]; 

% Distribuição de Coeficiente de Seção
  Cx = SAC./(2*W.*(T.-K));  % [1]     Coef. de Seção  
  
% Geração de Superfície
    N_points   = 5*ceil(0.5*B +T);  % [1] Pontos na baliza
    
    Cx_0      = 0.5; 
    Cx_switch = 0.6244;  
    Cx_1      = 1;
    
    Y_surf = zeros(N_points,length(x));
    Z_surf = zeros(N_points,length(x));
   
  for n = 1:length(x) 
    if (Cx(n) <= Cx_1)&&(Cx(n) > Cx_switch)
      Z_surf(:,n) = buntSectionIt(W(n),K(n),T,Cx(n),N_points)(:,1);
      Y_surf(:,n) = buntSectionIt(W(n),K(n),T,Cx(n),N_points)(:,2);
      
    elseif (Cx(n) <= Cx_switch)&&(Cx(n) >= Cx_0)
      Z_surf(:,n) = slimSectionIt(W(n),K(n),T,Cx(n),N_points)(:,1);
      Y_surf(:,n) = slimSectionIt(W(n),K(n),T,Cx(n),N_points)(:,2);
    endif  
  endfor
  
% Grandezas Derivadas
    rho = 1025;                             % [kg/m3] Densidade da água salgada
    
    dispV   = trapz(x,SAC);                 % [m3]    Deslocamento Volumétrico 
    dispM   = rho*dispV;                    % [kg]    Deslocamento Mássico
    Lcb     = trapz(x,x.*SAC)/dispV;        % [m]     Coord. Long. Centro de Carena (absoluta a partir da popa)
    Lcb_rel = (Lcb -Lbp/2)/Lbp;             % [1]     Coord. Long. Centro de Carena (relativa à meia nau)
    Vcb     = f_VCB(x,Y_surf,Z_surf);       % [m]     Coord. Vert. Centro de Carena (absoluta a partir da Linha de Base)
    Cb      = dispV/(Lbp*B*T);              % [1]     Coef. de Bloco 
    Cp      = Cb/Cm;                        % [1]     Coef. Prismático
    Awl     = trapz(x,2*W);                 % [m2]    Área do Plano de Linha D'água
    S       = surfArea(x,Y_surf,Z_surf);    % [m2]    Área Molhada
    Cwp     = Awl/(B*Lbp);                  % [1]     Coef. de Plano de Linha D'água 
    Lcf     = trapz(x,x.*(2*W))/Awl;        % [m]     Coord. Long. Baricentro do Plano de Linha D'água (absoluta a partir da popa)
    Lcf_rel = (Lcf -0.5*Lbp)/Lbp;           % [1]     Coord. Long. Baricentro do Plano de Linha D'água (relativa à meia nau)
    IzxWL   = trapz(x,(2*W).^3/12);         % [m4]    Momento de Inércia Transversal do Plano de Linha D'água

% Conferência
    format short
    printf('dispV   = %s m^3 \n',       num2str(dispV))
    printf('dispM   = %s ton \n',       num2str(1e-3*dispM))
    printf('Lcb     = %s m \n',         num2str(Lcb))
    printf('Lcb_rel = %s percent of Lbp from amidships \n',...
                                        num2str(100*Lcb_rel))
    printf('Vcb     = %s m \n',         num2str(Vcb))
    printf('Cb      = %s \n',           num2str(Cb))
    printf('Cp      = %s \n',           num2str(Cp))
    printf('Cwp     = %s \n',           num2str(Cwp))
    printf('Awl     = %s m^2 \n',       num2str(Awl))
    printf('IzxWL   = %s m^4 \n',       num2str(IzxWL))
    printf('S       = %s m^2 \n',       num2str(S))
    printf('Cx_min  = %s \n',           num2str(min(Cx)))
    printf('Cx_max  = %s \n',           num2str(max(Cx)))
  
  figure 1
    subplot(4,1,1)
      plot(x,SAC,'linewidth',1.5,'b',...
           zeros(length(x),1) + X_aft, linspace(0,1.2*max(SAC),length(x)),'-.k',...
           zeros(length(x),1) + X_fwd, linspace(0,1.2*max(SAC),length(x)),'-.k')
      hold on
      grid on
      xlim([0 Lbp])
      ylim([-0.1 1.2 ]*max(SAC))  
      ylabel('SAC [m2]')
      xlabel('x [m]')
      
    subplot(4,1,2)
      plot(x,K,'linewidth',1.5,'b',...
           x,zeros(length(x),1),':k','linewidth',1,...
           x,zeros(length(x),1) .+T,':k','linewidth',1,...
           zeros(length(x),1) + X_aft, linspace(-2.5,1.15*T,length(x)),'-.k',...
           zeros(length(x),1) + X_fwd, linspace(-2.5,1.15*T,length(x)),'-.k')
      text(min([2 0.05*Lbp]), max([-1 -0.25*T]),'BL','fontsize',7)
      text(max([Lbp-3 0.9*Lbp]),min([1.5*T T+1]), 'DWL','fontsize',7)
      grid on
      hold on
      xlim([0 Lbp]) 
      ylim([-1 (1.2*T+2)]) 
      xlabel('x [m]')
      ylabel('K [m]')
    
    subplot(4,1,3)
      plot(x,W,'linewidth',1.5,'b',...
           x,-W,':b','linewidth',0.75,...
           x,zeros(length(W),1),':k',...
           zeros(length(x),1) + X_aft, linspace(-0.6*B,0.6*B,length(x)),'-.k',...
           zeros(length(x),1) + X_fwd, linspace(-0.6*B,0.6*B,length(x)),'-.k')
           text(1, min([2 0.15*B]),'CL','fontsize',7)
      grid on
      hold on
      xlim([0 Lbp])
      ylim([-0.65*B 0.65*B])
      ylabel('W [m]')
      xlabel('x [m]')
      
    subplot(4,1,4)
      plot(x,Cx,'linewidth',1.5,'b',...
           x,zeros(length(x),1) +1,':r','linewidth',1,...
           x,zeros(length(x),1) +0.5,':r','linewidth',1,...
           zeros(length(x),1) + X_aft, linspace(0.1,1.1*max(Cx),length(x)),'-.k',...
           zeros(length(x),1) + X_fwd, linspace(0.1,1.1*max(Cx),length(x)),'-.k')
      grid on
      xlim([0 Lbp])
      ylim([0 1.2*max(Cx)])
      ylabel('Cx [1]')
      xlabel('x [m]')  
  
% Seleção de Balizas para plotagem  
    aft_ratio = (X_aft)/Lbp; 
    pmb_ratio = (X_fwd -X_aft)/Lbp;
    fwd_ratio = (Lbp -X_fwd)/Lbp;
    
    N_sections_aft = floor(1*aft_ratio*length(x));
    N_sections_pmb = floor(1*pmb_ratio*length(x));
    N_sections_fwd = floor(1*fwd_ratio*length(x));
    N_sections     = N_sections_aft ...   % [1] Qtde de balizas para desenho
                    +N_sections_pmb ...
                    +N_sections_fwd;     
    
    X_sections0_aft = linspace(0,X_aft,N_sections_aft);
    X_sections0_pmb = linspace(X_aft+1,X_fwd-1,N_sections_pmb);
    X_sections0_fwd = linspace(X_fwd,Lbp,N_sections_fwd);
    
    X_sections0 = horzcat(X_sections0_aft,X_sections0_pmb,X_sections0_fwd);   % [m] Coord. Long. em cada baliza
    Y_sections  = zeros(N_points,N_sections);                                 % [m] Coord. Transv. em cada baliza
    Z_sections  = zeros(N_points,N_sections);                                 % [m] Coord. Vert. em cada baliza
  
  for m = 1:N_sections
    x_idx = find(x == round(X_sections0(m)));
    X_sections(m)   = x(x_idx);
    Y_sections(:,m) = Y_surf(:,x_idx);
    Z_sections(:,m) = Z_surf(:,x_idx);  
  endfor

  figure 2
    for r = 1:N_sections
      plot3(zeros(N_points,1) +X_sections(r),Y_sections(:,r),Z_sections(:,r),'linewidth',2.5,'k',...
            zeros(N_points,1) +X_sections(r),-Y_sections(:,r),Z_sections(:,r),'color',[0,0,0]+0.5)
      hold on
    endfor
    
    plot3(x,W,zeros(length(x),1) +T,'linewidth',2.5,'b',...
          x,-W,zeros(length(x),1) +T,'linewidth',1.5,':b',...
          x,zeros(length(x),1),K,'linewidth',2.5,'r')
    xlim([0 Lbp])
    ylim([-1.1*max(W),1.1*max(W)])
    zlim([0,1.1*T])
    axis("equal")
    grid on

  Output = zeros(N_points,4, N_sections);

    for n = 1:N_sections
    Output(:,1,n) = zeros(N_points,1) +n;
    Output(:,2,n) = zeros(N_points,1) +X_sections(n);
    Output(:,3,n) = Y_sections(:,n);
    Output(:,4,n) = Z_sections(:,n);
  endfor
 
  fid = fopen('hullGeometry.txt', 'w+');
  
  fprintf(fid, '%f', 0);
  fprintf(fid, '\n');

  for n = 1:size(Output, 3)
    fprintf(fid, '\n');
    for m = 1:size(Output, 1)
      fprintf(fid, '%f ', Output(m,:,n));
      fprintf(fid, '\n');
    endfor
  endfor
  
  fclose(fid);
  
##  for n = 1:N_sections
##    Output(:,1,n) = zeros(N_points,1) +n+1;
##    Output(:,2,n) = zeros(N_points,1) +X_sections(n);
##    Output(:,3,n) = Y_sections(:,n);
##    Output(:,4,n) = Z_sections(:,n);
##  endfor
##  
##  Station0 = [(zeros(N_points,1))+1,...
##              (zeros(N_points,1) -0.01),...
##              (zeros(N_points,1)),...
##              (linspace(K(1),T,N_points))'];
##  
##  Station1 = [(zeros(N_points,1)+(N_sections +2)),...
##              (zeros(N_points,1)+Lbp+10),...
##              (zeros(N_points,1)),...
##              (linspace(K(end),T,N_points))'];
##  
##  fid = fopen('hullGeometry.txt', 'w+');
##  
##  fprintf(fid, '%f', 0);
##  fprintf(fid, '\n');
##  
##  for m = 1:size(Station0,1)
##    fprintf(fid, '\n');
##    fprintf(fid, '%f ', Station0(m,:))
##  endfor 
##
##  fprintf(fid, '\n');
##  
##  for n = 1:size(Output, 3)
##    fprintf(fid, '\n');
##    for m = 1:size(Output, 1)
##      fprintf(fid, '%f ', Output(m,:,n));
##      fprintf(fid, '\n');
##    endfor
##  endfor
##  
##  for m = 1:size(Station1,1)
##    fprintf(fid, '\n');
##    fprintf(fid, '%f ', Station1(m,:));
##  endfor 
##  
##  fclose(fid);
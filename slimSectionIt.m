function slimSecIt = slimSectionIt(W,K,T,Cx,N_points)
    sk0 = 0;
    sk1 = 3;
    E = 1;
    count = 0;
    
    while (E > 1E-9)&&(count < 100)
      sk_mid = mean([sk0,sk1]);
      section_mid = slimSection(sk_mid,N_points);
      z_mid = section_mid(:,1);
      y_mid = section_mid(:,2);
      Cx_mid = trapz(z_mid,y_mid);
      
      if Cx_mid == Cx
        break
        
      elseif  Cx_mid > Cx
        sk1 = sk_mid;
        
      elseif  Cx_mid < Cx
        sk0 = sk_mid;  
      endif
      
      E = abs((Cx_mid -Cx)/Cx);
      count += 1;
    endwhile
    
    sk = sk_mid;
    z0 = slimSection(sk,N_points)(:,1); 
    y0 = slimSection(sk,N_points)(:,2); 
    
    slimSecIt(:,1) = (T -K)*z0 +K;
    slimSecIt(:,2) = W*y0;
  endfunction
  
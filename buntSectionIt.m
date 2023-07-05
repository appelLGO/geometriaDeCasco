function buntSecIt = buntSectionIt(W,K,T,Cx,N_points)
    s0 = 1;
    s1 = 200;
    E = 1;
    count = 0;
    
    while (E > 1E-9)&&(count < 100)
      s_mid = mean([s0,s1]);
      section_mid = buntSection(s_mid,N_points);
      z_mid = section_mid(:,1);
      y_mid = section_mid(:,2);
      Cx_mid = trapz(z_mid,y_mid);
      
      if Cx_mid == Cx
        break
        
      elseif  Cx_mid > Cx
        s1 = s_mid;
        
      elseif  Cx_mid < Cx
        s0 = s_mid;  
      endif
      
      E = abs((Cx_mid -Cx)/Cx);
      count += 1;
    endwhile
    
    s = s_mid;
    z0 = buntSection(s,N_points)(:,1); 
    y0 = buntSection(s,N_points)(:,2);
    
    buntSecIt(:,1) = (T -K)*z0 +K;
    buntSecIt(:,2) = W*y0;
  endfunction

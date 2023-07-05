function A = surfArea(x,y,z)
  M = length(x);
  N = size(y)(1);
  
  
  a = zeros(M,1);
  
  for m = 1:(M -1)
    for n = 1:(N -1)
      P1 = [x(m),   y(n,m),     z(n,m)]';
      P2 = [x(m+1), y(n,m+1),   z(n,m+1)]';
      P3 = [x(m),   y(n+1,m),   z(n+1,m)]';
      P4 = [x(m+1), y(n+1,m+1), z(n+1,m+1)]';
      
      Pl = P2 -P1;
      Pt = P3 -P1;
      Pd = P4 -P1;
      
      S1 = cross(Pl,Pd)./2;
      S2 = cross(Pd,Pt)./2;
      A1 = sqrt(dot(S1,S1));
      A2 = sqrt(dot(S2,S2));
      
      a(m) += A1 +A2;
    endfor
  endfor  
  A = 2*sum(a);
endfunction
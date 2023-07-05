function VCB = f_VCB(x,y,z)
   for m = 1:length(x)
    A(m)  = trapz(z(:,m),2*y(:,m));
    mz(m) = trapz(z(:,m),2*z(:,m).*y(:,m));
  endfor
    Vol = trapz(x,A);
    VCB = trapz(x,mz)/Vol;
endfunction
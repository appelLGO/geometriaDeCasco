function buntSec = buntSection(s,N_points)
    z = linspace(0,1,N_points);
    y = ((1 .-(1-z).^s).^(1/s));
    
    buntSec(:,1) = z;
    buntSec(:,2) = y;
  endfunction
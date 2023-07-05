function slimSec = slimSection(slope_keel,N_points)
    slope_slim_keel = slope_keel;     % [m/m]  Deriv. Casco na Quilha       
    z = linspace(0,1, N_points);
   
    y = spline([0, 1],[slope_slim_keel, 0, 1, 0], z);
   
    slimSec(:,1) = z;
    slimSec(:,2) = y;   
  endfunction
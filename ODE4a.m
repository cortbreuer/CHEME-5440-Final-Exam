function dxdt = ODE4a(t, x, kd)
 
dxdt(1, 1) = 2 - (x(1) * x(2)) - (kd * x(1));

dxdt(2, 1) = (x(1) * x(2)) - x(2);
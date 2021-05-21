function dxdt = ODE3a(t, x)
 
dxdt(1, 1) = x(2) - (x(1) * x(2));

dxdt(2, 1) = (x(1) * x(2)) - x(2);


 

 



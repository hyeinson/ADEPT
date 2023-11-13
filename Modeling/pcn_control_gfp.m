function [t, y] = pcn_control_gfp(tspan, y0)

options = odeset('NonNegative',[1,2,3]);

[t,y] = ode45(@eom, tspan, y0, options);

end


function dydt = eom(Time,y)
    
global u_s alpha delta eta u_p gamma d


dydt = [ u_s * y(1) * (1 - (y(1)+y(2))) + alpha * y(2) - delta  * y(1) - eta * y(1) * y(2);
        
         u_p * y(2) * (1 - (y(1)+y(2))) - alpha * y(2) + eta * y(1) * y(2) ;
        
         gamma * y(2) - d * y(3);
       ];
   
  
end

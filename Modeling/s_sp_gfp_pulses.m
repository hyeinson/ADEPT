

function [t, y] = s_sp_gfp_pulses(tspan, y0, overall_time, alpha_duration)

    options = odeset('NonNegative', [1, 2, 3]);
    
    [t, y] = ode45(@(t, y) eom(t, y, overall_time, alpha_duration), tspan, y0, options);  % Pass overall_time to eom
end

function dydt = eom(Time, y, overall_time, alpha_duration)
    global u_s u_p alpha eta delta gamma d 
   
    internal_time = overall_time + Time; 
    
    if internal_time >= (12 - alpha_duration) && internal_time < 12
        cutting = 1;
    elseif internal_time >= (36 - alpha_duration) && internal_time < 36
        cutting = 1;
    elseif internal_time >= (60 - alpha_duration) && internal_time < 60
        cutting = 1;
    else 
        cutting = 0;
    end
    
    if internal_time > 12 && internal_time < 24
        killing = 1;
    elseif internal_time > 36 && internal_time < 48
        killing = 1;
    elseif internal_time > 60 && internal_time < 72
        killing = 1;
    else
        killing = 0;
    end
    
    
    dydt = [u_s * y(1) * (1.0 - (y(1) + y(2))) + alpha * y(2) * cutting - delta * y(1) * killing - eta * y(1) * y(2);
            u_p * y(2) * (1.0 - (y(1) + y(2))) - alpha * y(2) * cutting + eta * y(1) * y(2);
            gamma * y(2) - d * y(3);
           ];
end


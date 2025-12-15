function total_power_per_location=total_wind_power()
    rho=1.225;
    areas=[30,50,70];
    Cps=[0.35,0.4,0.45];
    rated_powers=[15000,30000,50000];
    wind_speeds=[
        5.0, 6.0, 5.5, 7.0, 6.5; % Location 1
        6.0, 7.5, 6.5, 8.0, 7.0; % Location 2
        4.0, 4.5, 5.0, 5.5, 6.0  % Location 3
    ];
    total_power_per_location=zeros(1, 3);
    for loc=1:3
        location_total_power=0;
        for turb=1:3
            for hour=1:5
                v=wind_speeds(loc, hour);
                cal_pwr=0.5*rho*areas(turb)*Cps(turb)*v^3;
                if cal_pwr>rated_powers(turb)
                    power_for_hour=rated_powers(turb);
                else
                    power_for_hour=cal_pwr;
                end
                location_total_power=location_total_power+power_for_hour;
            end
        end
        total_power_per_location(loc)=location_total_power;
    end
end
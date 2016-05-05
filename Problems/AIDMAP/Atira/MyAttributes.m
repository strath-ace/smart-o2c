classdef MyAttributes
    % The atributes of each node should be noted in this file
    
    properties  
        dV_dep %dV for departure
        dV_arr %dV needed at arrival
        dV_tot %total dV
        kep_trans %The keplerian elements describing the s/c orbit (the transfer)
        r_dep %departure r
        v_dep %departure velocity
        r_arr %arrival r
        v_arr %arrival v
        t_dep %departure time
        t_arr %arrival time
        tof %time of flight
        tof_tot %total tof so far
        body_dep %departure body
        body_arr %arrival body
        lambertV_ini %initial lambert velocity
        lambertV_final %final lambert velocity
    end

                
    
end
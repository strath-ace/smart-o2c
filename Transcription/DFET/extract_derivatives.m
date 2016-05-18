function [dx,du,dxb] = extract_derivatives(J,structure)

nx = (structure.state_order+1)*structure.num_eqs*structure.num_elems;               % number of state coefficients
nu = (structure.control_order+1)*structure.num_controls*structure.num_elems;        % number of control coefficients
nb = sum(structure.free_final_states);                                              % number of free final states (for DFET)

dx = zeros(size(J,1),nx);
du = zeros(size(J,1),nu);
dxb = zeros(size(J,1),nb);

startmat = 1;
startx = 1;
startu = 1*(nb>0);
startxb = 1*(nb>0);

for i = 1:structure.num_elems

    endx = startx+(structure.state_order+1)*structure.num_eqs-1;
    endmat = startmat+(structure.state_order+1)*structure.num_eqs-1;
    dx(:,startx:endx) = J(:,startmat:endmat);
    startx = endx+1;

    startmat = endmat+1*(structure.num_controls>0);
    endmat = startmat+(structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    
    endu = startu + (structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    % ugly, but apparently no empty to empty assignment can be done
    
    if nu>0
    
        du(:,startu:endu) = J(:,startmat:endmat);
        startu = endu+1;
        
    end
    
    
    startmat = endmat+1;
    
end

if startxb>0
   
    dxb = J(:,startmat:end);
    
else
    
    dxb = [];
    
end

end
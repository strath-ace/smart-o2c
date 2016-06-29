function [c,ceq,Jc,Jceq] = Simple_dynamics (f,structure,x_in,x_0,x_f,t_0,dfx,dfu)

t_f = x_in(1);
els = structure.uniform_els;
%els = x_in(2:2+structure.num_elems-2);
    
%els = [[0; els] [els;1]];
%x = x_in(2+structure.num_elems-1:end);

x = x_in(2:end);

c = [];
Jc = [];

[ceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,els,0,dfx,dfu);

% [ceq,Jceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,els,1,dfx,dfu);
% 
% %ceq = ceq.*structure.scale_transcribed_vars(2:end)/structure.scale_transcribed_vars(1);
% 
% Jceq = Jceq';
% 
% f2 = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f+1e-6,els,0,[],[]);
% dfdt = (f2-ceq)'/1e-6;
% 
% Jextra = dfdt;
% 
% % for i = 1:structure.num_elems-1
% %    
% %     models = els;
% %     models(i,2) = models(i,2)+1e-6;
% %     models(i+1,1) = models(i+1)+1e-6;
% %     f2 = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,models,0,[],[]);
% %     dfdtau = (f2-ceq)'/1e-6;
% %     Jextra = [Jextra;dfdtau];
% %     
% % end
% 
% Jceq = [Jextra;Jceq];
% 
end
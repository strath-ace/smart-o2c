close all
clear
clc

% load memory file

load('memory3.mat');

system('del dtlz.txt');                  % win

mem = memory;


hypervols = zeros(size(mem));
for i = 1:size(mem,1)   % problem (ZDT1,2,3,4)
    
    for j = 1:size(mem,2)   % num objectives (3,6,8)
        
        hv_str = [];
        
        no = 3*(j==1)+6*(j==2)+8*(j==3);
        nd = 5*(i==1)+10*(i>1)+no-1;
        
        hv_ref = ones(1,no)*1.5*(i==1)+ones(1,no)*2*(i>1);
        
        % write file
        
        for k = 1:size(mem,3)   % instance
            
            % extract points of the front
            qq = mem{i,j,k};
            points = qq(:,nd+1:nd+no);
            max(points);
            
            points_clean = points(all(points<hv_ref,2),:);
            
            if size(points_clean,1)<size(points,1)
                
                points = points_clean;
                
            end
            
            
            % write file to compute HV
            dlmwrite('dtlz.txt','#','-append');
            dlmwrite('dtlz.txt',points,'delimiter',' ','precision',9,'-append');
            
        end
        
        dlmwrite('dtlz.txt','#','-append');
        
        %         if j==1
        %           figure()
        %           plot3(points(:,1),points(:,2),points(:,3),'b.')
        %           size(points);
        %         end
        %
        % create executable string
        wfg_exec = 'wfg';%'/home/la/Dropbox/smart-o2c/Misc/WFG_1.15/wfg';
        curfile = 'dtlz.txt';%'/home/la/Dropbox/smart-o2c/Problems/MACS/DTLZ/dtlz.txt';
        
        for k = 1:length(hv_ref)-1
            
            hv_str = [hv_str num2str(hv_ref(i)) ', '];
            
        end
        
        hv_str = [hv_str num2str(hv_ref(end))];
        
        out_redir = '>HV.txt';
        
        str_to_call = [wfg_exec ' ' curfile ' ' hv_str ' ' out_redir];
        
        % call HV external code to compute it
        system(str_to_call);
        
        % remove temporary dtlz.txt file
        %system('rm dtlz.txt');     % linux
        system('del dtlz.txt');                  % win
        
        % extract HV from file
        hypervols(i,j,:) = load('HV.txt')./hv_ref(1)^no;
        hv_str = [];
        
    end
    
end

hv_mean = mean(hypervols,3);
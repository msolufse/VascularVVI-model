function Run_DriverBasic
clear all
close all
clc

for i = 1:2
        if i == 1
            call_flag  = 'PH'; 
        elseif i == 2
            call_flag  = 'Control'; 
        end
        
        switch call_flag
            case 'PH';
                [pars,~,~,~,data] = IPAH;
            case 'Control';
                [pars,~,~,~,data] = CONTROL;
        end

        for j = 1:3
            if j == 1
                data.model   = 'Nonlinear';
                data.INDMAP  = [1 2 7 8 9 10 11 13 17 18];  
                data.Fig     = 1;
            elseif j == 2
                data.model   = 'Linear';
                data.INDMAP  = [1 2 7 8 9 10 11 12 13 14]; 
                data.Fig     = 2;
            else
                data.model   = 'NoVVI';
                data.INDMAP  = [1 2 7 8 9 10 11 12 13 14];
                data.Fig = 3;
            end
            
            task = 1;     % Run LM Opt
            %task = 2;    % Plot Model
            
            par = 0; % Nom pars org (0), LM Opt pars (2)
            
            R = 3;   % Residual vector from manuscript
            
            data.R = R;
            
            if par == 2
                s = strcat('LMOpt_',call_flag,'_',data.model,'.mat');
                load(s);

                pars(history.INDMAP) = optpars(history.INDMAP);   % LM returns opt pars in array with all pars
            end
            
            if task == 1 % Optimization
                [history, optpars] = DriverBasic_opt(data,pars);
                s = strcat('Opt_',call_flag,'_',data.model,'.mat');
                save(s,'history','optpars','pars','data');
                
            elseif task == 2 % Plotting Model
                [rout,J,CV] = DriverBasic_plot(pars,data,call_flag,i,jj);
             
            end

        end
    end

end

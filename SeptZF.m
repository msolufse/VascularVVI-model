
function Vspt0 = SeptZF(Vlv,Vrv,Vdspt,pars)
    % Elatance Heart parameters
    EMrv   = pars(11); % max left ventricle elatance
    Emrv   = pars(12); % min left ventricle elatance
    EMlv   = pars(13); % max left ventricle elatance
    Emlv   = pars(14); % min left ventricle elatance

    % Heart timing parameters
    Tcrv = pars(15);        % right ventricle contraction
    Trrv = Tcrv + pars(16); % right ventricle relaxation
    Tclv = Tcrv;            % left ventricle contraction MITCHEL SCALE!
    Trlv = Trrv;  
    
    Llv   = pars(17);
    Lrv   = pars(18);
    P0rv  = pars(22);
    P0lv  = pars(23);

    % Septum Parameters
    EMspt  = pars(19);  % Septum FW elstnce (mmHg/mL)
    Emspt  = pars(20);  % Septum ED pressure param (mmHg)
    Lspt   = pars(21);  % Septum ED pressure param (1/mL)
    P0spt  = pars(24);
    
    % Heart elastance
    erv  = ElastanceDriver(0,Tcrv,Trrv);   % Right ventricle elastance
    elv  = ElastanceDriver(0,Tclv,Trlv);   % Right ventricle elastance
    espt = ElastanceDriver(0,Tclv,Trlv);
    
    V0 = 5;
        
    % Calculate the zero function which is derived from the fact that
    %  the pressure attributed to the septal wall is the difference between
    %  the pressure attributed to the left ventricular free wall minus the 
    %  pressure attributed to the left ventricular free wall --> 
    %  P_spt = P_lvf - P_rvf
    
    % We will use the complete linear model and analytical Vspt solution to
    % find VSpt0 the initial condition

    Vspt0  = ((elv'*EMlv + (1-elv')*Emlv).*(Vlv - V0) - (erv'*EMrv + (1-erv')*Emrv).*(Vrv - V0) ...
            + (espt'*EMspt + (1-espt')*Emspt)*Vdspt)./(espt'*EMspt + (1-espt')*Emspt ...
            + elv'*EMlv + (1-elv')*Emlv + erv'*EMrv + (1-erv')*Emrv);

end


function [rout,J,CV] = CVmodelLn(pars,data)


global ODE_TOL REL_TOL parsUN

VsaU   = parsUN(1);
VsvU   = parsUN(2);
VpaU   = parsUN(3);
VpvU   = parsUN(4);
VvU    = parsUN(5);
VsptU  = parsUN(6);

pars = exp(pars);

Init = data.Init;

% Resistances
Rs    = pars(1); % systemic periphery
Rp    = pars(2); % pulmonary periphery
Rava  = pars(3); % aortic valve
Rmva  = pars(4); % mitral valve
Rpva  = pars(5); % pulmonary valve
Rtva  = pars(6); % tricuspid valve

% Compliances
Csa = pars(7);  % systemic arteries
Csv = pars(8);  % systemic veins
Cpa = pars(9);  % pulmonary arteries
Cpv = pars(10); % pulmonary veins

% Elastance Heart parameters
EMrv   = pars(11); % max right ventricle elatance
Emrv   = pars(12); % min right ventricle elastance
EMlv   = pars(13); % max left ventricle elatance
Emlv   = pars(14); % min left ventricle elastance

% Heart timing parameters
Tcrv = pars(15);        % right ventricle contraction
Trrv = Tcrv + pars(16); % right ventricle relaxation

Tclv = Tcrv;            % left ventricle contraction
Trlv = Trrv;            % left ventricle relaxation

% VVI parameters
Llv   = pars(17);
Lrv   = pars(18);
P0lv  = pars(23);
P0rv  = pars(22);

% Septum Parameters
EMspt  = pars(19);  % Septum FW elstnce (mmHg/mL)
Emspt  = pars(20);  % Septum ED pressure param (mmHg)
Lspt   = pars(21);  % Septum ED pressure param (1/mL)
P0spt  = pars(24);

% Arteries and venous pressures
CV.ppaS   = [];  
CV.pPAdS  = [];
CV.ppaMS  = [];
CV.pPAMdS = [];    
CV.ppamS  = [];
CV.pPAmdS = [];

CV.psaS   = [];
CV.psaMS  = [];
CV.pSAMdS = [];
CV.psamS  = [];
CV.pSAmdS = [];

CV.ppvS   = [];
CV.ppvmS  = [];
CV.ppvdS  = [];
CV.psvS   = [];
CV.psvmS  = [];
CV.psvdS  = [];

% Heart pressures
CV.psptS  = [];
CV.plvS   = [];
CV.plvMS  = [];
CV.pLVMdS = [];
CV.plvmS  = [];
CV.pLVmdS = [];

CV.prvS   = [];
CV.prvMS  = [];
CV.pRVMdS = [];
CV.prvmS  = [];
CV.pRVmdS = [];

% Flow
CV.COS    = [];
CV.COdS   = [];

% Heart Volumes
CV.VlvS   = [];
CV.VlvMS  = [];
CV.VlvMdS = [];
CV.VlvmS  = [];
CV.VlvmdS = [];
CV.VrvS   = [];
CV.VrvMS  = [];
CV.VrvMdS = [];
CV.VrvmS  = [];
CV.VrvmdS = [];
CV.VsptS  = [];

CV.VlvfS  = [];
CV.VrvfS  = [];

% Arterial and venous volumes
CV.VsaS   = [];
CV.VsameS = [];
CV.VsadS  = [];
CV.VsvS   = [];
CV.VsvmeS = [];
CV.VsvdS  = [];
CV.VpaS   = [];
CV.VpameS = [];
CV.VpadS  = [];
CV.VpvS   = [];
CV.VpvmeS = [];
CV.VpvdS  = [];

% Time
CV.tdS    = [];
CV.tCOS   = [];

T  = data.T;
td = data.td;  
dt = data.dt; 
NC = 1; 

k1 = 1; % index of first time step in first period
k2 = round(T/dt)+k1; %index of last time step in first period
notconv = true;
notdone = true;
tdS     = data.td;
NCend = 0;
while notdone && NC < 100  
    clear Vsa Vsv Vpa Vpv Vrv Vlv Vspt
    clear psa psv ppa ppv Erv prv Elv plv pspt
    clear qs qp  qmva qtva 
        
    tdc = data.td + (NC-1)*T;
%   disp([tdc(1) tdc(end)]);
    
    options=odeset('RelTol',REL_TOL, 'AbsTol',ODE_TOL);   
    Init = Init(1:6);

    sol = ode15s(@modelBasic_VVILn,tdc,Init,options,pars,tdc(1),data); %solves the ODE uses modelBasic
    sols= deval(sol,tdc);   

    % extract solutions
    Vsa  = sols(1,:)';  % systemic arteries
    Vsv  = sols(2,:)';  % systemic veins
    Vpa  = sols(3,:)';  % pulmonary arteries
    Vpv  = sols(4,:)';  % pulmonary veins
    Vrv  = sols(5,:)';  % right ventricle
    Vlv  = sols(6,:)';  % left ventricle
        
    psa = (Vsa-VsaU)/Csa; % systemic arteries
    psv = (Vsv-VsvU)/Csv; % systemic veins
    ppa = (Vpa-VpaU)/Cpa; % pressure pulmonary artery
    ppv = (Vpv-VpvU)/Cpv; % pressure pulmonary vein

    % Heart elastance
    for w = 1:length(tdc) 
        elv(w)  = ElastanceDriver(tdc(w)-tdc(1),Tclv,Trlv); % Left ventricle elastance
        erv(w)  = ElastanceDriver(tdc(w)-tdc(1),Tcrv,Trrv); % Right ventricle elastance
        espt(w) = ElastanceDriver(tdc(w)-tdc(1),Tclv,Trlv);
    end
   
    Vspt  = ((elv'*EMlv  + (1-elv')*Emlv).*(Vlv - VvU) - (erv'*EMrv + (1-erv')*Emrv).*(Vrv - VvU) ...
          + (espt'*EMspt + (1-espt')*Emspt)*VsptU)./(espt'*EMspt + (1-espt')*Emspt ...
          + elv'*EMlv + (1-elv')*Emlv + erv'*EMrv + (1-erv')*Emrv);

    % Ventricles
    Vrvf = (Vrv - VvU) + Vspt;
    Vlvf = (Vlv - VvU) - Vspt;

    prv = erv'.*EMrv.*(Vrvf) + (1-erv').*Emrv.*(Vrvf);
    plv = elv'.*EMlv.*(Vlvf) + (1-elv').*Emlv.*(Vlvf);
    
    % Linear Septum
    pspt = espt'.*EMspt.*(Vspt - VsptU) + (1-espt').*Emspt.*(Vspt - VsptU);
     
    % Flows defined by Ohm's Law
    qs  = (psa - psv)/Rs; % systemic periphery
    
    % Cardiac Output
    COs  = trapz(tdc, qs)/(tdc(end)-tdc(1))/1000*60;
    
    % time
    CV.tdS  = [CV.tdS tdc(1:end-1)]; 
    
    % pspt
    CV.psptS   = [CV.psptS pspt(1:end-1)']; 
    
    % plv
    CV.plvS   = [CV.plvS   plv(1:end-1)'];
    CV.plvMS  = [CV.plvMS  max(plv)];
    CV.pLVMdS = [CV.pLVMdS data.pSlv];
    CV.plvmS  = [CV.plvmS  min(plv)];
    CV.pLVmdS = [CV.pLVmdS data.pDlv];

    % prv
    CV.prvS   = [CV.prvS   prv(1:end-1)'];
    CV.prvMS  = [CV.prvMS  max(prv)];
    CV.pRVMdS = [CV.pRVMdS data.pSrv];
    CV.prvmS  = [CV.prvmS  min(prv)];
    CV.pRVmdS = [CV.pRVmdS data.pDrv];
    
    % psv
    CV.psvS   = [CV.psvS   psv(1:end-1)'];    
    CV.psvmS  = [CV.psvmS  mean(psv)];
    CV.psvdS  = [CV.psvdS  data.psv];
    
    % ppa
    CV.ppaS   = [CV.ppaS   ppa(1:end-1)'];    
    CV.ppaMS  = [CV.ppaMS  max(ppa)];  
    CV.pPAMdS = [CV.pPAMdS data.pPAM];    
    CV.ppamS  = [CV.ppamS  min(ppa)];
    CV.pPAmdS = [CV.pPAmdS data.pPAm];
    
    % psa
    CV.psaS   = [CV.psaS   psa(1:end-1)'];
    CV.psaMS  = [CV.psaMS  max(psa)];
    CV.pSAMdS = [CV.pSAMdS data.pSAM];
    CV.psamS  = [CV.psamS  min(psa)];
    CV.pSAmdS = [CV.pSAmdS data.pSAm];

    % ppv
    CV.ppvS   = [CV.ppvS    ppv(1:end-1)'];
    CV.ppvmS  = [CV.ppvmS   mean(ppv)];
    CV.ppvdS  = [CV.ppvdS   data.ppv];
    
    % CO
    CV.COS    = [CV.COS COs];
    CV.COdS   = [CV.COdS data.CO]; 
    
    % Volumes
    CV.VsaS   = [CV.VsaS   Vsa(1:end-1)'];
    CV.VsameS = [CV.VsameS mean(Vsa)];
    CV.VsadS  = [CV.VsadS  data.Vsa];
    
    CV.VsvS   = [CV.VsvS   Vsv(1:end-1)'];
    CV.VsvmeS = [CV.VsvmeS mean(Vsv)];
    CV.VsvdS  = [CV.VsvdS  data.Vsv];
    
    CV.VpaS   = [CV.VpaS   Vpa(1:end-1)'];
    CV.VpameS = [CV.VpameS mean(Vpa)];
    CV.VpadS  = [CV.VpadS  data.Vpa];
    
    CV.VpvS   = [CV.VpvS   Vpv(1:end-1)'];
    CV.VpvmeS = [CV.VpvmeS mean(Vpv)];
    CV.VpvdS  = [CV.VpvdS  data.Vpv];
    
    % Heart Volumes
    CV.VlvS   = [CV.VlvS   Vlv(1:end-1)'];
    CV.VlvMS  = [CV.VlvMS  max(Vlv)];
    CV.VlvMdS = [CV.VlvMdS data.VlvD];
    CV.VlvmS  = [CV.VlvmS  min(Vlv)];
    CV.VlvmdS = [CV.VlvmdS data.VlvS];

    CV.VrvS   = [CV.VrvS   Vrv(1:end-1)'];
    CV.VrvMS  = [CV.VrvMS  max(Vrv)];
    CV.VrvMdS = [CV.VrvMdS data.VrvD];
    CV.VrvmS  = [CV.VrvmS  min(Vrv)];
    CV.VrvmdS = [CV.VrvmdS data.VrvS];
     
    CV.VlvfS  = [CV.VlvfS   Vlvf(1:end-1)'];
    CV.VrvfS  = [CV.VrvfS   Vrvf(1:end-1)'];    
    CV.VsptS  = [CV.VsptS   Vspt(1:end-1)'];
    
    Init = [Vsa(end) Vsv(end) Vpa(end) Vpv(end) Vrv(end) Vlv(end) Vspt(end)];
    
    if NC > 20
        dPsaM  = CV.psaMS(end) - CV.psaMS(end-1);
        dPsam  = CV.psamS(end) - CV.psamS(end-1);
        dPpaM  = CV.ppaMS(end) - CV.ppaMS(end-1);
        dPpam  = CV.ppamS(end) - CV.ppamS(end-1);
        dPrvM  = CV.prvMS(end) - CV.prvMS(end-1);
        dPrvm  = CV.prvmS(end) - CV.prvmS(end-1);
        
        diff  = max([dPsaM dPsam dPpaM dPpam dPrvM dPrvm]);
        
        if diff < 1e-4
          notconv = false;
        end
        
        if notconv
            NCend = 0; 
            % to make sure that if for whatever reason we arent in steady 
            % state again we restart
        else
            NCend = NCend + 1;
        end
        
        if NCend > 5
           notdone = false;
        end   
    end
    
    NC = NC + 1;
end

% save solution (end point)
CV.tdS    = [CV.tdS    tdc(end)];
CV.psptS  = [CV.psptS  pspt(end)];
CV.plvS   = [CV.plvS   plv(end)];
CV.prvS   = [CV.prvS   prv(end)];
CV.psvS   = [CV.psvS   psv(end)];
CV.ppaS   = [CV.ppaS   ppa(end)];
CV.psaS   = [CV.psaS   psa(end)];
CV.ppvS   = [CV.ppvS   ppv(end)];

CV.VsaS   = [CV.VsaS   Vsa(end)];
CV.VsvS   = [CV.VsvS   Vsv(end)];
CV.VpaS   = [CV.VpaS   Vpa(end)];
CV.VpvS   = [CV.VpvS   Vpv(end)];
CV.VlvS   = [CV.VlvS   Vlv(end)];
CV.VrvS   = [CV.VrvS   Vrv(end)];
CV.VsptS  = [CV.VsptS  Vspt(end)];

CV.VlvfS  = [CV.VlvfS   Vlvf(end)];
CV.VrvfS  = [CV.VrvfS   Vrvf(end)];


M2   = round(NC-5);

routS = [(CV.plvMS(M2:end) -CV.pLVMdS(M2:end))./CV.pSAMdS(M2:end)/sqrt(M2) ...
         (CV.psaMS(M2:end) -CV.pSAMdS(M2:end))./CV.pSAMdS(M2:end)/sqrt(M2) ... 
         (CV.psamS(M2:end) -CV.pSAmdS(M2:end))./CV.pSAMdS(M2:end)/sqrt(M2) ...
         (CV.prvMS(M2:end) -CV.pRVMdS(M2:end))./CV.pPAMdS(M2:end)/sqrt(M2) ...
         (CV.ppaMS(M2:end) -CV.pPAMdS(M2:end))./CV.pPAMdS(M2:end)/sqrt(M2) ...
         (CV.ppamS(M2:end) -CV.pPAmdS(M2:end))./CV.pPAMdS(M2:end)/sqrt(M2) ...
 	     (CV.prvmS(M2:end) -CV.pRVmdS(M2:end))./CV.ppvdS(M2:end)/sqrt(M2) ... 
         (CV.psvmS(M2:end) -CV.psvdS(M2:end)) ./CV.ppvdS(M2:end)/sqrt(M2) ...
         (CV.plvmS(M2:end) -CV.pLVmdS(M2:end))./CV.ppvdS(M2:end)/sqrt(M2) ... 
         (CV.ppvmS(M2:end) -CV.ppvdS(M2:end)) ./CV.ppvdS(M2:end)/sqrt(M2) ...
         2*(CV.COS(M2:end)   -CV.COdS(M2:end))  ./CV.COdS(M2:end) /sqrt(M2)  ...
         (CV.VlvMS(M2:end) -CV.VlvMdS(M2:end))./CV.VlvMdS(M2:end)/sqrt(M2) ...
         (CV.VlvmS(M2:end) -CV.VlvmdS(M2:end))./CV.VlvMdS(M2:end)/sqrt(M2) ...
         (CV.VrvMS(M2:end) -CV.VrvMdS(M2:end))./CV.VlvMdS(M2:end)/sqrt(M2) ...
         (CV.VrvmS(M2:end) -CV.VrvmdS(M2:end))./CV.VlvMdS(M2:end)/sqrt(M2) ]; 

rout = routS';
J    = rout'*rout;
CV.M2 = NC-1;
W_c = 1.33322*10^(-4);
CV.work = [trapz(plv,Vlv) trapz(prv,Vrv)].*W_c;
CV.EF   = (max([Vlv Vrv])-min([Vlv Vrv]))./max([Vlv Vrv]);


end


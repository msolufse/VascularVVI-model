function xdot = modelBasic_VVILn(t,y,pars,tst,data)

global parsUN

Vsa    = y(1);  % systemic artery volume
Vsv    = y(2);  % systemic venous volume
Vpa    = y(3);  % pulmonary artery volume
Vpv    = y(4);  % pulmonary venous volume
Vrv    = y(5);  % right ventricular volume
Vlv    = y(6);  % left ventricular volume

% Unstressed Volumes
VsaU  = parsUN(1);
VsvU  = parsUN(2);
VpaU  = parsUN(3);
VpvU  = parsUN(4);
VvU   = parsUN(5);
VsptU = parsUN(6);

% Resistances
Rs   = pars(1);  % systemic arteries
Rp   = pars(2);  % pulmonary arteries
Rava = pars(3);  % aortic valve 
Rmva = pars(4);  % mitral valve
Rpva = pars(5);  % pulmonary valve
Rtva = pars(6);  % tricuspid valve


% Compliances
Csa = pars(7);  % systemic arteries
Csv	= pars(8);  % systemic veins
Cpa = pars(9);  % pulmonary arteries
Cpv = pars(10); % pulmonary veins

% Heart
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

% Septum Parameters
EMspt  = pars(19);  % Septum FW elstnce (mmHg/mL)
Emspt  = pars(20);  % Septum ED pressure param (mmHg)
Lspt   = pars(21);  % Septum ED pressure param (1/mL)
P0spt  = pars(24);

% Timevarying elastance
elv   = ElastanceDriver(t-tst,Tclv,Trlv);  % Elastance left ventricle
erv   = ElastanceDriver(t-tst,Tcrv,Trrv);  % Elastance right ventricle
espt  = ElastanceDriver(t-tst,Tclv, Trlv); 

% Pressure
psa  = (Vsa-VsaU)/Csa; % systemic arteries 
psv  = (Vsv-VsvU)/Csv; % systemic veins 
ppa  = (Vpa-VpaU)/Cpa; % pressure pulmonary artery
ppv  = (Vpv-VpvU)/Cpv; % pressure pulmonary vein

Vspt  = ((elv'*EMlv  + (1-elv')*Emlv).*(Vlv - VvU) - (erv'*EMrv + (1-erv')*Emrv).*(Vrv - VvU) ...
      + (espt'*EMspt + (1-espt')*Emspt)*VsptU)./(espt'*EMspt + (1-espt')*Emspt ...
      + elv'*EMlv + (1-elv')*Emlv + erv'*EMrv + (1-erv')*Emrv);

% Ventricles
Vrvf = (Vrv - VvU) + Vspt;
Vlvf = (Vlv - VvU) - Vspt;

prv = erv'.*EMrv.*(Vrvf) + (1-erv').*Emrv.*(Vrvf);
plv = elv'.*EMlv.*(Vlvf) + (1-elv').*Emlv.*(Vlvf);

% Valves
if plv > psa % flow through aortic valve 
   qava = (plv-psa)/Rava; 
else
   qava  = 0;
end

if prv > ppa % flow through pulmonary valve
   qpva = (prv-ppa)/Rpva;
else
   qpva  = 0;
end

% Peripheral flows
qs   = (psa - psv) / Rs; % Systemic periphery
qp   = (ppa - ppv) / Rp; % Pulmonary periphery

% Venous/Valve flows (with venous valves - flow cannot go backwards into the veins)
if psv > prv  % flow through tricuspid valve
    qtva = (psv-prv)/Rtva; % valve open
else
    qtva = 0;             % valve closed
end

if ppv > plv % flow through mitral valve
    qmva = (ppv-plv)/Rmva; % valve open
else
    qmva = 0;             % valve closed
end

% Differential Equations
dVsa = qava - qs;
dVsv = qs   - qtva;
dVpa = qpva - qp;
dVpv = qp   - qmva;
dVrv = qtva - qpva;
dVlv = qmva - qava;

xdot = [dVsa;...
        dVsv;...
        dVpa;...
        dVpv;...
        dVrv;...
        dVlv];
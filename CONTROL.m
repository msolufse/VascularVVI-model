% This function initializes the parameters for the model and sets initial
% values for the variables.
function [pars, Init,low,hi,data] = CONTROL

global DIFF_INC ODE_TOL REL_TOL

ODE_TOL  = 1e-10;
REL_TOL  = 1e-10;
DIFF_INC = sqrt(ODE_TOL);
data.DIFF_INC = DIFF_INC;
data.ODE_TOL  = ODE_TOL;
data.REL_TOL  = REL_TOL;

% Define subject specific quantities
BSA = 1.85; % Made same as PH patient

G = 2;
if G == 2
  TotalVol = (3.47*BSA - 1.954)*1000; % Female
elseif G == 1
  TotalVol = (3.29*BSA - 1.229)*1000; % Male
end

% Patient-specific Cardiac Output
CO   = 5;              % Cardiac output (L/min)
qtot = CO*1000/60;     % Flow (mL/sec)

% Pressures (related to subject) 
% Systemic arteries (mmHg)
pSsa  = 120;  % systolic systemic artery pressure               % 
pDsa  = 80;   % diastolic systemic artery pressure 
psa   = 2/3*pDsa+1/3*pSsa;  % mean systemic artery pressure

% Systemic Veins (mmHg)
psv    = 4;  % mean systemic vein pressure 

% Pulmonary arteries (mmHg)
pSpa  = 21;   % systolic pulmonary artery pressure
ppa   = 12;   % mean pulmonary arterial pressure
pDpa  = (ppa-1/3*pSpa)*3/2;  % diastolic pulmonary arterial pressure

% Pulmonary Veins (mmHg)
ppv    = 6;  % mean pulmonary vein pressure 

% Ventricles (mmHg)
pMlv   = 1.025*pSsa;  % maximum left ventricular pressure
pmlv   = 0.975*ppv;   % minimum left ventricular pressure
pMrv   = 1.025*pSpa;  % maximum right ventricular pressure
pmrv   = 0.975*psv;   % minimum right ventricular pressure

% Heart Rate (b/min)
Ave_HR = 65;

% Timing parameters 
dt     = 0.01;   % Time step seconds (resolution saved within each cycle)
T      = round(60/Ave_HR/dt)*dt; % Length of cardiac cycle (seconds)
NC     = 1;      % Run code for at least one cardiac cycle
td     = 0:dt:T; % Number of points per cardiac cycle

Tcrv = 0.30; % ventricular contraction
Trrv = 0.20; % ventricular relaxation

% Ventricular Volumes Data Bredfelt (mL)
SV   = CO/Ave_HR*1000;  % Stroke volume 

VMrv = 120;          % Max right ventricular volume 
VMlv = 120;          % Max left ventricular volume
Vmrv = VMrv - SV;    % min right ventricular volume                    
Vmlv = VMlv - SV;    % min left ventricular volume
percentRV = VMrv/TotalVol; % Percent of total volume in right ventricle
percentLV = VMlv/TotalVol; % Percent of total volume in left ventricle
                  
% Volume distribution (Boron Text) (mL)
Vsa  = 0.13*TotalVol; % Arterial volume (13%)
Vsv  = (0.715-percentRV-percentLV)*TotalVol;% Venous volume (65%)
Vpa  = 0.03*TotalVol; % Pulmonary artery (3%)
Vpv  = 0.125*TotalVol; % Pulmonary vein volume (11%)

% Unstressed Volumes (from Beneken) (mL)
VsaU = Vsa*(1-0.27);  % Stressed volume 27% of arterial volume
VsvU = Vsv*(1-0.075); % Stressed volume 7.5% of venous volume
VpaU = Vpa*(1-0.58);  % Stressed volume 58% of pulmonary artery 
VpvU = Vpv*(1-0.11);  % Streseed volume 11% of pulmonary vein 
VvU  = 10;            % Unstressed ventricular volume

% Resistences (mmHg s/mL)
% Peripheral resistances (Ohm's law)
Rs   = (psa-psv)/qtot;   % systemic resistance
Rp   = (ppa-ppv)/qtot;   % pulmonary resistance

% Valve resistances (mmHg s/mL)
Rava = (pMlv-pSsa)/qtot; % resistance aortic valve
Rpva = (pMrv-pSpa)/qtot; % resistance pulmonary valve
Rmva = (ppv-pmlv)/qtot;  % resistance mitral valve
Rtva = (psv-pmrv)/qtot;  % resistance tricuspid valve 

% Compliances (mL/mmHg)
Csa = (Vsa-VsaU)/pSsa;    % systemic artery compliance
Csv = (Vsv-VsvU)/psv;     % systemic venous compliance  
Cpa = (Vpa-VpaU)/pSpa;    % pulmonary artery compliance
Cpv = (Vpv-VpvU)/ppv;     % pulmonary vein compliance

% Heart parameters (Elastance - E (mmHg/mL))
% Elastance parameters ventricles
Emrv = pmrv/(VMrv-VvU); % min right ventricle
EMrv = pMrv/(Vmrv-VvU); % max right ventricle 
Emlv = pmlv/(VMlv-VvU); % min left ventricle
EMlv = pMlv/(Vmlv-VvU); % max left ventricle
EMmf = Emlv/EMlv;

% Nonlinear model (lambda's and P0s) ventricles
Lrvf = 0.007; % lambda right ventricle
Llvf = 0.007; % lambda left ventricle

P0rv = Emrv/Lrvf;  % P0 right ventricle set to agree with linear model
P0lv = Emlv/Llvf;  % P0 left ventricle set to agree with linear model

% VVI Parameters   
% Septum free wall parameters 
% VVI Parameters   
% Septum free wall parameters 
EMspt  = 10*EMlv;          % Septum FW elstnce (mmHg/mL) 
Vdspt  = 0;                % Septum zero P volume (mL)
Emspt  = EMmf*EMspt;       % Septum ES elastance
Lspt   = 50*Llvf;          % Septum ED elastance 
P0spt  = Emspt/Lspt;       % Septum zero pressure related to create agreement with linearized model

% Parameter vector
x0 = [Rs Rp Rava Rmva Rpva Rtva  ...       % 1-6
      Csa Csv Cpa Cpv ...                  % 7-10
      EMrv Emrv EMlv Emlv ...              % 11-14
      Tcrv Trrv Llvf Lrvf ...              % 15-18         
      EMspt Emspt Lspt P0rv P0lv P0spt]';  % 19-24        
      
global parsUN
parsUN = [VsaU VsvU VpaU VpvU VvU Vdspt];      
data.parsUN = parsUN;

% Bounds for optimization
% hi  = x0*1.25; %pars + log(4);%
% low = x0*0.75; %pars - log(4);%
% 
% hi  = log(hi);
% low = log(low);

pars = log(x0);
hi  = pars + log(4);
low = pars - log(4);

% Timing parameters upper and lower bounds (not estimated)
hi(15:16,:)  = pars(15:16,:) + log(2);
low(15:16,:) = pars(15:16,:) - log(2);


% Initial septal volume
Vspt = SeptZF(VMlv,VMrv,Vdspt,x0);

Vrvf = VMrv - VvU;
Vlvf = VMlv - VvU;

% Initial conditions
Init = [Vsa Vsv Vpa Vpv Vrvf Vlvf Vspt]; 

% set up data structure

% Timing parameters
data.T      = T;
data.NC     = NC;
data.dt     = dt;
data.td     = td;

% Patient characteristics
data.Gender = G;
data.BSA    = BSA;
data.HR     = Ave_HR;
data.CO     = CO; 

% Pressure (data)
data.pSAM   = pSsa;
data.pSAm   = pDsa;
data.pPAM   = pSpa;
data.pPAm   = pDpa;

data.ppv    = ppv;
data.psv    = psv;
data.pSrv   = pMrv;
data.pDrv   = pmrv;
data.pSlv   = pMlv;
data.pDlv   = pmlv;

% Volumes
data.VlvS   = Vmlv;
data.VlvD   = VMlv;
data.VrvS   = Vmrv;
data.VrvD   = VMrv;
data.Vsa    = Vsa;
data.Vsv    = Vsv;
data.Vpa    = Vpa;
data.Vpv    = Vpv;
  
% Initial conditions and bounds
data.Init = Init;
data.hi   = hi;
data.low  = low;
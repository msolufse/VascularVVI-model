function [history, x] = DriverBasic_opt(data,pars)
%par_cluster = 0;

gPars.ALLPARS = pars;   % all parameters
gPars.INDMAP  = data.INDMAP; 

optx   = pars(gPars.INDMAP);

%set up structure for optimization
data.ALLPARS  = pars;

opthi     = data.hi(data.INDMAP);
optlow    = data.low(data.INDMAP);

% LM Kelley
% optimization information
maxiter = 50;        % max number of iterations        
mode    = 2;         % Performs Levenberg-Marquart optimization
nu0     = 2.d-1;     % Regularization parameter
[xopt, histout, costdata, jachist, xhist, rout, sc] = ...
    newlsq_v2(optx,'opt_wrap',1.d-3,maxiter,mode,nu0,...
              opthi,optlow,data);
history.x      = xhist;
history.hist   = histout;

history.INDMAP   = data.INDMAP;

x = data.ALLPARS;
x(data.INDMAP) = xopt;


end



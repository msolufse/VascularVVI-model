%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function [rout,J,CVsave] = model_wrap(pars,data)

tpars = data.ALLPARS;
tpars(data.INDMAP') = pars;

model = data.model;

switch model
    case 'Nonlinear'
        [rout,J,CVsave] = CVmodelNL(tpars,data);
    case 'Linear'
        [rout,J,CVsave] = CVmodelLn(tpars,data);
    case 'NoVVI'
        [rout,J,CVsave] = CVmodelNoVVI(tpars,data);
end

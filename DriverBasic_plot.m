function [rout,J,CV] = DriverBasic_plot(pars,data,call_flag,I)

model = data.model;
FS    = 14;
switch model
    case 'Nonlinear'
        [rout,J,CV] = CVmodelNL(pars,data);
        mod = 1;
        my_col = repmat(linspace(0,0,10)',1,3); 
        line = 2;
        name = 'Nonlinear';
    case 'Linear'
        [rout,J,CV] = CVmodelLn(pars,data);
        mod = 2;
        my_col = repmat(linspace(1,1,10)',1,3); 
        my_col(:,1) = 0;
        line = 2;
        name = 'Linear';
    case 'NoVVI'
        [rout,J,CV] = CVmodelNoVVI(pars,data);
        mod = 3;
        my_col = repmat(linspace(1,1,10)',1,3); 
        my_col(:,2:3) = 0;
        line = 2;
        name = 'No VVI';
end

N    = length(data.td);
tdc  = data.td;

plv  = CV.plvS(end-N+1:end);
psa  = CV.psaS(end-N+1:end);
psv  = CV.psvS(end-N+1:end);
prv  = CV.prvS(end-N+1:end);
ppa  = CV.ppaS(end-N+1:end);
ppv  = CV.ppvS(end-N+1:end);
CO   = CV.COS(end); 
try
    pspt = CV.psptS(end-N+1:end);
    Vspt = CV.VsptS(end-N+1:end);
catch
end
Vrv  = CV.VrvS(end-N+1:end);
Vlv  = CV.VlvS(end-N+1:end);
Vpa  = CV.VpaS(end-N+1:end);
Vpv  = CV.VpvS(end-N+1:end);
Vsa  = CV.VsaS(end-N+1:end);
Vsv  = CV.VsvS(end-N+1:end);

s = strcat(call_flag);


hFig = figure(I);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [200 200 800 800]);

 sgtitle(s,'fontsize',FS+4)

maxpLV = data.pSlv*ones(size(tdc));
minpLV = data.pDlv*ones(size(tdc));
subplot(4,4,1); hold on;
  plot(tdc,plv,'color', my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  plot(tdc,maxpLV,'--k','Linewidth',line,'HandleVisibility','off');
  plot(tdc,minpLV,'--k','Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Plv');
  ylim([0 140]);

maxpRV = data.pSrv*ones(size(tdc));
minpRV = data.pDrv*ones(size(tdc));
subplot(4,4,2); hold on
  plot(tdc,prv,'color', my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  plot(tdc,maxpRV,'--k','Linewidth',line,'HandleVisibility','off');
  plot(tdc,minpRV,'--k','Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Prv');
  ylim([0 data.pSrv*1.1]);

subplot(4,4,3); hold on
  try
      plot(tdc,pspt,'color', my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  catch
  end
  set(gca,'FontSize',FS);
  ylabel('Pspt');
  ylim([0 110]);

maxvLV = data.VlvD*ones(size(tdc));
minvLV = data.VlvS*ones(size(tdc));
subplot(4,4,5); hold on
  plot(tdc,Vlv,'color', my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  plot(tdc,maxvLV,'--k','Linewidth',line,'HandleVisibility','off');
  plot(tdc,minvLV,'--k','Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Vlv');
  ylim([40 130]);

maxvRV = data.VrvD*ones(size(tdc));
minvRV = data.VrvS*ones(size(tdc));
subplot(4,4,6); hold on
  plot(tdc,Vrv,'color', my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  plot(tdc,maxvRV,'--k','Linewidth',line,'HandleVisibility','off');
  plot(tdc,minvRV,'--k','Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Vrv');
  ylim([data.VrvS*0.95 data.VrvD*1.05]);

subplot(4,4,7); hold on
  try
    plot(tdc,Vspt,'color', my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  catch
  end
  set(gca,'FontSize',FS);
  ylabel('Vspt');
  ylim([-1 10]);

subplot(4,4,13); hold on
  plot(Vlv,plv,'color', my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',16);
  ylabel('Plv');
  xlabel('Vlv');
  ylim([0 140]);
  xlim([40 130]);

subplot(4,4,14); hold on
  plot(Vrv,prv,'color',my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Prv');
  xlabel('Vrv');
  xlim([data.VrvS*0.95 data.VrvD*1.05]);
  ylim([0 data.pSrv*1.1]);

subplot(4,4,15); hold on
try
  plot(Vspt,pspt,'color',my_col(1,:),'Linewidth',line,'HandleVisibility','off');
catch
end
set(gca,'FontSize',FS);
xlabel('Vspt');
ylabel('Pspt');
ylim([0 110]);
xlim([-1 10]);

maxpSA = data.pSAM*ones(size(tdc));
minpSA = data.pSAm*ones(size(tdc));
subplot(4,4,9); hold on
   plot(tdc,psa,'color',my_col(1,:),'Linewidth',line,'HandleVisibility','off');
   plot(tdc,maxpSA,'--k','Linewidth',line,'HandleVisibility','off');
   plot(tdc,minpSA,'--k','Linewidth',line,'HandleVisibility','off');
   set(gca,'FontSize',FS);
   ylabel('psa');
   xlabel('time');
   ylim([75 125]);

meanpsv = data.psv*ones(size(tdc));
subplot(4,4,10); hold on
  plot(tdc,psv,'color',my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  plot(tdc,meanpsv,'--k','Linewidth',2,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Psv');
  xlabel('time');
  ylim([0.9*data.psv 1.1*max(psv)])
  
maxpPA = data.pPAM*ones(size(tdc));
minpPA = data.pPAm*ones(size(tdc));
subplot(4,4,11);hold on;
  plot(tdc,ppa,'color',my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  plot(tdc,maxpPA,'--k','Linewidth',line,'HandleVisibility','off');
  plot(tdc,minpPA,'--k','Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Ppa');
  xlabel('time');
  ylim([data.pPAm*0.9 data.pPAM*1.1]);

minppv = data.ppv*ones(size(tdc)); 
subplot(4,4,12);hold on
  plot(tdc,ppv,'color',my_col(1,:),'Linewidth',line,'HandleVisibility','off');
  plot(tdc,minppv,'--k','Linewidth',line,'HandleVisibility','off');
  set(gca,'FontSize',FS);
  ylabel('Ppv');
  xlabel('time');
  ylim([0.9*data.ppv 1.1*max(ppv)])

COt  = CO*ones(size(tdc));
COdt = data.CO*ones(size(tdc));
subplot(4,4,16);hold on;
  plot(tdc,COdt,'--k','Linewidth',line,'HandleVisibility','off');
  plot(tdc,COt,'color',my_col(1,:),'Linewidth',line,'DisplayName',name);
  ylim([min(data.CO-0.5,min(COt)) max(data.CO+0.5,max(COt))]);
  set(gca,'FontSize',FS);
  ylabel('CO');
  xlabel('time');

ax = subplot(4,4,4,'Visible','off');
axPos = ax.Position;
delete(ax);

Lgnd = legend('show');
Lgnd.Position(1) = axPos(1);
Lgnd.Position(2) = axPos(2);

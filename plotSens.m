function plotSens(history,call,I)
%% Comment this out if excluding timing parameters
% 
% Rsens = history.Rsens;
% Isens = history.Isens;
% sens_norm = history.Nsens;
% sens  = history.sens;

Names = {'Rs','Rp','Rava','Rmva','Rpva','Rtva', ...
         'Csa','Csv','Cpa','Cpv',...
         'EMRv','EmRv','EMLv','EmLv',...
         'Tcrv','Trrv','Llvf','Lrvf',...
         'EMspt','Emspt','Lspt', ...
         'P0rv','P0lv','P0spt'};

sens  = history.sens;

% ranked classical sensitivities without timing
[M,N] = size(sens);
k = 1;
for i = [1:N]
    sens_norm(k)=norm(sens(:,i),2);
    k = k+1;
end

[Rsens,Isens] = sort(sens_norm,'descend');


Isens

figure ;clf;%hold on;
h=semilogy(sens_norm(Isens)./sens_norm(Isens(1)),'*');
set(h,'linewidth',2);
set(h,'Markersize',10);
set(gca,'Fontsize',16);
grid on;
ylabel('Sensitivities');
xlabel('Parameters')
xticks(1:length(Isens));
xlim([0,length(Isens)]);
xticklabels(Names(Isens));
xtickangle(45);
set(gca,'yscale','log')
print -dpng Sensitivities.png

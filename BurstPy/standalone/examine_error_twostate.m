%% BENCHMARK TO TEST INFERENCE METHODS
    % fit: 2 state models
    % varying params: sumk=kon+koff
    
T_range=[600 600];         % Trace length in seconds
ncell_range=[50 200];      % Number of cells
pon_range=0.1;                          % Probability of Pon
isrounded=0;
mkdir('Results');

for cnt=1:numel(T_range)
    T=T_range(cnt);
    ncell=ncell_range(cnt);
        if ~exist(['Results/twostate_' num2str(pon_range) '_T' num2str(T) '_ncell' num2str(ncell) '.mat'],'file')
            inferencetwostate(T,ncell,pon_range,isrounded);
        end
end
%% Check the results
% Plot the normal condition
pon_range=0.1;
T=600;
ncell=50;
load(['Results/twostate_' num2str(pon_range) '_T' num2str(T) '_ncell' num2str(ncell) '.mat'],'sumk_range','sumkhat','itrmax');
figure;
h1=plot(sort(sumk_range),sort(sumk_range),'k--');hold on;
for i=1:numel(sumkhat)
    sumkhat{i}=sumkhat{i}(sumkhat{i}<0.15);
end
    sumkmean=cellfun(@mean,sumkhat);
    sumkerr=sqrt(cellfun(@var,sumkhat));
    %h2=errorbar(sumk_range+0.001,sumkmean,sumkerr);hold on;
    h2=shadedErrorBar(sumk_range+0.001,sumkmean,sumkerr,'b',0.1);hold on;
    set(h2.mainLine,'LineWidth',2);
%% Plot the ideal condition
T=600;
ncell=200;
load(['Results/twostate_' num2str(pon_range) '_T' num2str(T) '_ncell' num2str(ncell) '.mat'],'sumk_range','sumkhat','itrmax');
for i=1:numel(sumkhat)
    sumkhat{i}=sumkhat{i}(sumkhat{i}<0.15);
end
    sumkmean=cellfun(@mean,sumkhat);
    sumkerr=sqrt(cellfun(@var,sumkhat));
    %h3=errorbar(sumk_range,sumkmean,sumkerr);hold on;
    h3=shadedErrorBar(sumk_range+0.001,sumkmean,sumkerr,'r',0.1);hold on;
    set(h3.mainLine,'LineWidth',2);
xlabel('k_{on} + k_{off} (1/s)');
ylabel('Inferred k_{on} + k_{off} (1/s)');
axis([0 0.15 0 0.17]);
%%
legend([h2.mainLine h3.mainLine],'Condition1','Condition2');
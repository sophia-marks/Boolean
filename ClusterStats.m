%ClusterStats.m

clear;

load WolframClusterData.txt
wo = WolframClusterData;
load WuenscheClusterData.txt
wu = WuenscheClusterData;
[d1,d2] = size(wo); %same as size(wu)


plot(wo(:,2:d2),'LineWidth',2)
title('Wolfram vs. Wuensche', 'FontSize',12, 'FontWeight', 'bold')
hold on
plot(wu(:,2:d2),'--','LineWidth',2)
title('Wolfram vs. Wuensche cluster comparison', 'FontSize',12, 'FontWeight', 'bold')
legend('wo1', 'wo2', 'wo3', 'wo4', 'wu1', 'wu2', 'wu3', 'wu4','Location','best') %classes


%{
figure 
subplot(1,2,1)
boxplot(wo(:,2:d2))
title('Wolfram', 'fontsize', 12, 'FontWeight','bold');
subplot(1,2,2)
boxplot(wu(:,2:d2))
title('Wuensche', 'fontsize', 12, 'FontWeight','bold');
%}


[P1,ANOVATAB1,STATS1] = anova1(wo(:,2:d2))
figure
[COMPARISON1,MEANS1,H1,NAMES1] = multcompare(STATS1,'display','on');
num2cell(MEANS1)
COMPARISON1
%gtext(num2str(P1))

[P2,ANOVATAB2,STATS2] = anova1(wu(:,2:d2))
figure
[COMPARISON2,MEANS2,H2,NAMES2] = multcompare(STATS2,'display','on');
num2cell(MEANS2)
COMPARISON2
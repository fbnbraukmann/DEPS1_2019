%% piRNA analysis in deps-1
% load piRNA counts and total unique mapping read counts 
clear
load('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB100-FB199/FB169/dependent_libaries/FB169_pirna_boxplot.mat')

% label samples
samples = {'wt','wt','prg-1','prg-1','deps-1','deps-1'};

% normalise counts
normLibSizeCounts = bsxfun(@rdivide, piRNAcounts,totaluniquemappingreads)*10^6;

% calculate mean of the two biological samples
clear mean_count
for i = 1:2:5
    tmp = mean(normLibSizeCounts(:,i:i+1),2);
    mean_count(:,i) = tmp;
end

% remove empty columns and keep only means
mean_count(:,2) = [];
mean_count(:,3) = [];

% remove pirna for which no counts were observed
index1 = sum(mean_count,2) == 0;
mean_count(index1,:) = [];

% make a boxplot
boxplot(log(mean_count+1),'Symbol','k+','OutlierSize',4,'Colors','k')
strains = {'WT', 'prg-1', 'deps-1'};
% ylim([-0.5 10]);
xticks(1:length(strains));
xticklabels(strains);
ax = gca;
ax.LineWidth = 1.0;
ax.FontSize = 16
ax.FontName = 'Arial'
ax.Box = 'on'
ax.YLabel.String = 'log2(reads+1/million)';
ax.YLabel.FontSize = 16;
ax.YLabel.Position = [0.275 3.4215 -1.0000]
% line
% [h13, p13] = ttest(mean_count(:,1),mean_count(:,3));
% [h12, p12] =ttest(mean_count(:,1),mean_count(:,2));
% [h23, p23]  =ttest(mean_count(:,2),mean_count(:,3));
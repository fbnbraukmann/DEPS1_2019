clc
clear all
load('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB200-FB299/FB216/matlab/FB216_repeats.mat')

repeatsCountTable = readtable('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/sRNA/fea_repeats/fea_repeatmask.grouped.xlsx');
sum_star = readtable('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/sRNA/sum_star/star.summary.txt');
totalreads = sum_star.Var2;
geneid_repeats = repeatsCountTable.Geneid;

samples = repeatsCountTable(:,7:end).Properties.VariableNames;
strains = {'GR1799_1', 'GR1799_2', 'GR1799_3', 'SX1316_1', 'SX1316_2', 'SX1316_3', 'SX1636_1', 'SX1636_2', 'SX1636_3', 'SX1888_1', 'SX1888_2', 'SX1888_3', 'SX3249_1', 'SX3249_2', 'SX3249_3', 'SX3266_1', 'SX3266_2', 'SX3266_3', 'SX3396_1', 'SX3396_2', 'SX3396_3', 'SX3397_1', 'SX3397_2', 'SX3397_3'};
sample_subselection = [1:24];



groups = zeros(length(strains),length(samples));
for i =  1:length(strains)
    groups(i,:) = strncmp(samples, strains{i},length(strains{i}));
end
%
counts_repeats = repeatsCountTable{:,7:end};
normLibSizeCounts = bsxfun(@rdivide, counts_repeats,totalreads');
indx = nanmean(normLibSizeCounts,2) > 5*10^-5;
normLibSizeCounts = normLibSizeCounts(indx,:);
geneid = repeatsCountTable.Geneid;
geneid = geneid(indx);

average_norm_reads = [];
for i = 1:length(groups)
    average_norm_reads(:,i) = mean(normLibSizeCounts(:,logical(groups(i,:))),2);
end

foldchange_repeats = bsxfun(@rdivide, average_norm_reads,average_norm_reads(:,4));
%
map = brewermap(64,'RdBu');
%allsamples_biologicalreplicates_1
%sample_subselectionreduced  = [1:3:24]

%allsamples_biologicalreplicates_3
%sample_subselectionreduced  = [1:24];


%wtprgmutdepsdepsGW_biologicalrep_1
%sample_subselectionreduced  = [4,7,10,13,19]

%wtprgmutdepsdepsGW_biologicalrep_3
sample_subselectionreduced = [4:6,7:9,10:12,13:15,19:21];



data = foldchange_repeats(:,sample_subselectionreduced);
cg = clustergram(log(data), 'RowLabels', geneid,...
                             'ColumnLabels', strains(sample_subselectionreduced)',...
                             'RowPdist', 'euclidean',...
                             'ColumnPdist', 'euclidean',...
                             'Colormap',map,...
                             'DisplayRange',6);           
T = array2table(log10(data),...
    'VariableNames',strains(sample_subselectionreduced),...
    'RowNames',geneid);
writetable(T,'/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/matlab/output/repeats/FB324_read_repeats_heatmap20180802.xlsx',...
    'WriteVariableNames',1,...
    'WriteRowNames',1);
                         
fgr = plot(cg)
saveas(fgr,'/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/matlab/output/repeats/FB324_read_repeats_heatmap20180802.eps','epsc2');
                           

%%
filename = 'normalised_reads_exons.xlsx';
xlswrite(filename,norm_reads_exon,1)
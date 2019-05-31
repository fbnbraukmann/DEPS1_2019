%% loads exon count raw data and gene_enrichment_classes
clearvars -except exonCountTable
close all 
clc

exonCountTable = readtable('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/sRNA/fea_exon/fea_exon.transcriptid.xlsx');
sum_star = readtable('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/sRNA/sum_star/star.summary.txt');
totalreads = sum_star.Var2;
geneid = exonCountTable.Geneid;

samples = exonCountTable(:,7:end).Properties.VariableNames;

strains = {'GR1799_1', 'GR1799_2', 'GR1799_3', 'SX1316_1', 'SX1316_2', 'SX1316_3', 'SX1636_1', 'SX1636_2', 'SX1636_3', 'SX1888_1', 'SX1888_2', 'SX1888_3', 'SX3249_1', 'SX3249_2', 'SX3249_3', 'SX3266_1', 'SX3266_2', 'SX3266_3', 'SX3396_1', 'SX3396_2', 'SX3396_3', 'SX3397_1', 'SX3397_2', 'SX3397_3'};


groups = zeros(length(strains),length(samples));
for i =  1:length(strains)
    groups(i,:) = strncmp(samples, strains{i},length(strains{i}));
end
%
counts = exonCountTable{:,7:end};
normLibSizeCounts = bsxfun(@rdivide, counts,totalreads');
indx = nanmean(normLibSizeCounts,2) > 5*10^-5 ;
normLibSizeCounts = normLibSizeCounts(indx,:);
geneid = exonCountTable.Geneid;
geneid = geneid(indx);

average_norm_reads = [];
for i = 1:length(groups)
    average_norm_reads(:,i) = mean(normLibSizeCounts(:,logical(groups(i,:))),2);
end

foldchange_exon = bsxfun(@rdivide, average_norm_reads,average_norm_reads(:,4));


clear indx1 nz i pseudoRefSample ratios sizeFactors

load('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB200-FB299/FB216/gene_enrich/FB216_gene_enrichment_classes.mat');

classesnames = {'alg34', 'csr1', 'ergo1', 'germline', 'hrde1', 'pirna', 'soma', 'wago', 'geneid'}
classes = {alg34, csr1, ergo1, germline, hrde1, pirna, soma, wago, geneid};

clear overlapmatrix

for i = 1:length(classes)
    overlapmatrix{i} = ismember(cellstr(geneid),classes{i});
    matches(i) = sum(overlapmatrix{i});
end

clear notoverlapmatrix

for i = 1:length(classes)
    notoverlapmatrix{i} = ~ismember(classes{i},cellstr(geneid));
    notmatches(i) = sum(notoverlapmatrix{i});
end

clear notfound
for i = 1:length(classes)
    tmp = classes{i};
    notfound{i} = tmp(notoverlapmatrix{i});
end


clear tmp tmp2 tmp3 geneid_classes_tmp geneid_classes
for i = 1:length(classes)
    indexofclasses = overlapmatrix{i};
    %remove lines with nan for wt sample 1, most of the lines that are nan
    %in wt sample 1 have nan everywhere
    tmp2 = foldchange_exon(indexofclasses,:);
    geneid_classes_tmp = geneid(indexofclasses,:);
    tmp3 = tmp2(~isnan(tmp2(:,1)),:);
    geneid_classes{i} = geneid_classes_tmp(~isnan(tmp2(:,1)),:);
    reads_classes{i} = tmp3;
end

% 1 = alg34
% 2 = csr1
% 3 = ergo1
% 4 = germline
% 5 = hrde1
% 6 = pirna
% 7 = soma
% 8 = wago
% 9 = geneid
for i = 1:9
    close all hidden
    classofinterested = i;
    classesnames{i};

    % sample_subselectionreduced  = [4,7,10,13,19] for wt, prg-1, mut-16,
    % deps-1 and deps-1deltaGW
    sample_subselectionreduced  = [4:6,7:9,10:12,13:15,19:21];
    clear tmp tmp2 tmp3
    % the following line will plot only a selection of samples
    tmp = reads_classes{1,i};
    tmp2 = tmp(:,sample_subselectionreduced);
    %remove row nan in wild tyoe
    nonnanindex = ~isnan(tmp2(:,3));
    data =  tmp2(nonnanindex,:);
    %set 0  to small value above 0
    data(data == 0) = min(data(data > 0));
    map = brewermap(64,'RdBu');
    cg = clustergram(log(data),...
        'ColumnLabels', samples(sample_subselectionreduced)',...
        'RowPdist', 'euclidean',...
        'ColumnPdist', 'euclidean',...
        'Colormap',map,...
        'DisplayRange',6);
    addTitle(cg, classesnames{i})
         
    T = array2table(log(data),...
        'VariableNames',strains(sample_subselectionreduced),...
        'RowNames',geneid_classes{classofinterested});
    writetable(T,['/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/matlab/output/FB324_read_' classesnames{i} '_heatmap20180802.xlsx'],...
        'WriteVariableNames',1,...
        'WriteRowNames',1);
% To make the changes programmatically it is necessary to grab a handle
% to to the figure associated with the clustergram and then manipulate the figure properties.
% The following code demonstrates how to obtain the figure handle:
% First, make all handles visible. This is necessary because clustergram objects...
% are created with 'HandleVisibility' property set to 'off':
    close all hidden
    
ax = plot(cg)
fgr = gcf;
    set(0,'ShowHiddenHandles','on')
    %Then, get all handles under the root object:
    allhnds = get(0,'Children');
    %Find the handles that correspond to clustergram objects
    cgfigidx = strmatch('Clustergram',get(allhnds,'Tag'));
    cffighnd = allhnds(cgfigidx);
    %Make the non-visible handles hidden again
    set(0,'showhiddenHandles','off')
    %Finally, if there is more than one clustergram object, operate on the last one in the list.
    if length(cffighnd)>1
        warning('More than one clustergram handle found.  Using most recent clustergram')
        cffighnd = cffighnd(end);
    end
    
    c = colorbar(ax,'Location','manual','Position',[0.85 0.125 0.025 0.60]);
    c.Label.String = 'log fold change';
    saveas(ax,['/Users/braukmann/Dropbox (ericmiskalab)/Experiments/FB300-FB399/FB324/matlab/output/FB324_read_' classesnames{i} '_heatmap20180802.eps'],'epsc2');
end

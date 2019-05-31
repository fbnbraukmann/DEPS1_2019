% load bam
% load fasta file
% count reads at starting nucleotide along the fasta file 
% count how many are on the histone
% normalise reads
% plot histogram
% 

%%%%% use the load line to load data instead of calculating from bam file
%%%%% again

% load('/Users/braukmann/Dropboxericmiskalab/Experiments/FB169/piRNAsensorplot/FB169_piRNAsensorplot.mat')

%%%%%
%%%%%

clear
clc

input_folder = dir('/Users/braukmann/Dropbox (ericmiskalab)/Experiments/Archive/FB100-FB199/FB169/piRNAsensorplot/bam');
counts = zeros(2153-489,length(input_folder));
for i = 3:1:length(input_folder)
    if sum(input_folder(i).name(end-3:end) == '.bam') == 4
    bam_input = {[input_folder(i).folder, '/', input_folder(i).name]};
    reads = bamread(bam_input{1},'Pmex_5_gfp_his_58',[489 2153]);
    
    index_fwdread = [reads.Flag]' == 0;
    position = double([reads.Position])';
    clear read_length
    
    for j = 1:length(reads)
        read_length(j) = str2double(reads(j).CigarString(1:2));
    end
    read_length = read_length';
    
    position(~index_fwdread) = position(~index_fwdread) + read_length(~index_fwdread);
    position = position - 489;
    
    number_of_histon_reads = sum(position >= 1400-489 & position <= 1772-489);
    
    number_of_reads_at_each_nucleotide = histcounts(position,2153-489);
    normalised_reads_at_each_nucelotide = number_of_reads_at_each_nucleotide/number_of_histon_reads;
    counts(:,i) = normalised_reads_at_each_nucelotide';
    end
end

%% plot, give a number from the input_folder
samples = [24,26,32,34,38,40,76,78];
sample_name = {'wt', 'wt', 'mut-16', 'mut-16', 'prg-1', 'prg-1', 'deps-1', 'deps-1'};
for i = 1:8   
    subplot(4,2,i)
    bar(counts(:,samples(i)),'BarWidth',10,'FaceColor','k')
    title(sample_name(i));
    ylabel('reads/H2B reads')
    ylim([0 1.5])
    ax = gca;
    ax.LineWidth = 1.0;
    ax.FontSize = 16;
    ax.FontName = 'Arial';
    ax.Box = 'on';
    ax.YLabel.FontSize = 16;
end
%%
samples = [24,32,38,76];
sample_name = {'wt', 'mut-16', 'prg-1', 'deps-1'};
for i = 1:4   
    subplot(1,4,i)
    bar([1:900],counts([1:900],samples(i)),'BarWidth',10,'FaceColor',[51/255 204/255 51/255])
    hold on
    bar([901:1280],counts([901:1280],samples(i)),'BarWidth',10,'FaceColor',[204/255 204/255 204/255])
    hold on
    bar([1281:1359],counts([1281:1359],samples(i)),'BarWidth',10,'FaceColor',[51/255 51/255 153/255])
    hold on
    bar([1360:1664],counts([1360:1664],samples(i)),'BarWidth',10,'FaceColor',[0/255 153/255 255/255])
    title(sample_name(i));
    ylabel('reads/H2B reads')
    ylim([0 1.5])
    ax = gca;
    ax.LineWidth = 1.0;
    ax.FontSize = 16;
    ax.FontName = 'Arial';
    ax.Box = 'on';
    ax.YLabel.FontSize = 16;
end


%Script for analysing data from white noise experiments 2014. Written by
%Varun Saravanan December 2014. Modified by Varun Saravanan September 2018
%to include bootstrap sampling with sampling of birds and then iterations.

clear
close all
clc

%P is the structure in which all bird parameters and data will be stored.

% bird_names:
% 1. gr90wh83
% 2. gr93yw43
% 3. gr98gy55
% 4. pk54pu4
% 5. or57rd126
% 6. gr93yw43 - post 6-OHDA leaion
% 7. gr98gy55 - post 6-OHDA lesion
% 8. pk54pu4 - post 6-OHDA lesion
% 9. lb30rd4
% 10. lb32rd15
% 11. pk27pu87
% 12. bk26gy38
% 13. lb32rd15 - post 6-OHDA lesion
% 14. pk27pu87 - post 6-OHDA lesion
% that may be called for analysis. For reference only. Actual call is made
% with corresponding numbers.

%Prelesion birds: birds_to_use = [1 2 3 4 5 9 10 11 12];
% same_birds = [1 3 5 6];
% diff_birds = 2:9;
%Post 6OHDA lesion birds: birds_to_use = [6 7 8 13 14];
%same_birds = 2;
%diff_birds = birds;
%Post sham lesion: birds_to_use = [15 16 17 18];
%same_birds = [1 2 4];
%diff_birds = [2 3 4];
%Pre 6OHDA matched birds: birds_to_use = [2 3 4 10 11];
%dist_to_use = 1:9;
%Pre Sham matched birds: birds_to_use = [1 5 9 12];
%dist_to_use = 3:10;
%Lukas's single syl headphones expts:
% birds_to_use = [19 20 21 22 23 24];

wn = 0;

if wn
    birds_to_use = [1 2 3 4 5 9 10 11 12];
else
    birds_to_use = [19 20 21 22 23 24];
end

 %Can specify which birds in the set of all bird names are to be included in generating figures.
plot_non_bootstrapped = 1;
plot_baseline = 0; %Toggle variable to control whether or not to plot 3 days of baseline.
%Set these:
num_base_days = 3;
if birds_to_use(1)==19
    num_wn_days = 14;
else
    num_wn_days = 3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load all bird parameters into structure P.

counter = 1;
for i = birds_to_use
    P(counter) = load_bird_params_2018(i);
    counter = counter + 1;
end


num_days = num_base_days + num_wn_days;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating mean pitch change over groups - targeted syllables, same type
%syllables and different syllables. This uses the mean pitch shift of each
%syllable and averages over all syllables.

mean_targeted_pitch_shift = zeros(num_base_days + num_wn_days,1);
mean_same_pitch_shift = zeros(num_base_days + num_wn_days,1);
mean_diff_pitch_shift = zeros(num_base_days + num_wn_days,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To calculate the mean and std dev, we use resampling. Mean is defined as
%the mean of the sample of means each obtained by resampling first the
%birds and second the syllable iterations for each bird and condensing all
%the data. Error bars are 95% CI of the sample of means, namely 1.96*SD of
%sample of means.


same_birds = [];
diff_birds = [];
target_birds = [];
for j = 1:length(birds_to_use)
    if ~isempty(P(j).same_syls)
        same_birds = [same_birds; j];
    end
    if ~isempty(P(j).diff_syls)
        diff_birds = [diff_birds; j];
    end
    for i = 1:length(P(j).target_syl)
        target_birds = [target_birds; [j P(j).target_syl(i)]];
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First to calculate the mean, we use the actual data. This is more accurate
%than using the mean of the population of bootstrapped means though the 2
%values should be close to each other given enough resampling.

for i = 1:num_days
    temp2 = [];
    for j = 1:length(birds_to_use)
        if birds_to_use(1) == 19
            temp = horzcat(P(j).data{i,P(j).target_syl});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,P(j).target_syl});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_targeted_pitch_shift(i) = mean(temp2);
    
    temp2 = [];
    for j = 1:length(same_birds)
        if birds_to_use(1) == 19
            temp = horzcat(P(same_birds(j)).data{i,P(same_birds(j)).same_syls});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(same_birds(j)).data{i,P(same_birds(j)).same_syls});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_same_pitch_shift(i) = mean(temp2);
    
    temp2 = [];
    for j = 1:length(diff_birds)
        if birds_to_use(1) == 19
            temp = horzcat(P(diff_birds(j)).data{i,P(diff_birds(j)).diff_syls});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(diff_birds(j)).data{i,P(diff_birds(j)).diff_syls});
            temp2 = vertcat(temp2,temp);
        end
    end
    mean_diff_pitch_shift(i) = mean(temp2);
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We also need to write a section to compute bootstrap samples for stats.
%These will involve 10000 samples and may need some sort of matrix to be
%computed so that averages can be taken.

%So, in the headphones paper, we made this easier by computing a 4-D
%bootstrap matrix from which to draw samples - so when we needed to average
%over multiple days, we just made a new column and drew from that final
%column. The sacrifice made to create such a matrix however was to treat
%each syllable as equally likely regardless of number of occurrence. While
%this was justified in the headphones paper, here the number of times each
%syllable is repeated matters...
%A way around this would be to create the same matrix as before but change
%the number of resamples drawn based on the number of iterations of that
%syllable. That way, we get the randomization desired and account for
%syllable occurrence.

%Bootstrapping part
nboot = 300; %No of times to resample for bootstrapping.
max_syls = 16;

bootstrapping_matrix = zeros(nboot,(num_days+1),max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            temp = P(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            bootstrapping_matrix(:,k,j,i) = datasample(temp,nboot);
        end
        for n = 1:nboot
            if birds_to_use(1) == 19 %Average over last 3 days for headphones
                bootstrapping_matrix(n,end,j,i) = nanmean(bootstrapping_matrix(n,(end-2):(end-1),j,i));
            else %Average over last two days for WN
                bootstrapping_matrix(n,end,j,i) = nanmean(bootstrapping_matrix(n,(end-2):(end-1),j,i));
            end
            temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
            if birds_to_use(1) == 19
                bootstrapping_matrix(n,:,j,i) = -1*P(i).shift_direction*12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            else
                bootstrapping_matrix(n,:,j,i) = P(i).shift_direction*12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
        end
    end
end

%Now do the bootstrapping to get error bars for figures: 
nboot2 = 1000; %Only repeating 1000 times for error bars
bootstats1 = zeros(nboot2,num_days);
bootstats2 = zeros(nboot2,num_days);
bootstats3 = zeros(nboot2,num_days);
for n = 1:nboot2
    
    temp_targ_birds = datasample(1:numel(birds_to_use),numel(birds_to_use));
    temp_data = [];
    for t = 1:length(temp_targ_birds)
        temp_syls = datasample(P(temp_targ_birds(t)).target_syl,length(P(temp_targ_birds(t)).target_syl));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,1:(end-1),temp_syls(s),temp_targ_birds(t)));
        end
    end
    bootstats1(n,:) = nanmean(temp_data,1);
    
    temp_same_birds = datasample(same_birds,length(same_birds));
    temp_data = [];
    for t = 1:length(temp_same_birds)
        temp_syls = datasample(P(temp_same_birds(t)).same_syls,length(P(temp_same_birds(t)).same_syls));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,1:(end-1),temp_syls(s),temp_same_birds(t)));
        end
    end
    bootstats2(n,:) = nanmean(temp_data,1);
    
    temp_diff_birds = datasample(diff_birds,length(diff_birds));
    temp_data = [];
    for t = 1:length(temp_diff_birds)
        temp_syls = datasample(P(temp_diff_birds(t)).diff_syls,length(P(temp_diff_birds(t)).diff_syls));
        for s = 1:length(temp_syls)
            temp_pulls = datasample(1:nboot,nboot);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,1:(end-1),temp_syls(s),temp_diff_birds(t)));
        end
    end
    bootstats3(n,:) = nanmean(temp_data,1);
end

std_targeted_pitch_shift = nanstd(bootstats1,0,1);
std_same_pitch_shift = nanstd(bootstats2,0,1);
std_diff_pitch_shift = nanstd(bootstats3,0,1);

if plot_baseline == 0
    mean_same_pitch_shift = mean_same_pitch_shift(3:end);
    mean_same_pitch_shift(1) = 0;
    std_same_pitch_shift = std_same_pitch_shift(3:end);
    std_same_pitch_shift(1) = 0;
    mean_diff_pitch_shift = mean_diff_pitch_shift(3:end);
    mean_diff_pitch_shift(1) = 0;
    std_diff_pitch_shift = std_diff_pitch_shift(3:end);
    std_diff_pitch_shift(1) = 0;
    mean_targeted_pitch_shift = mean_targeted_pitch_shift(3:end);
    mean_targeted_pitch_shift(1) = 0;
    std_targeted_pitch_shift = std_targeted_pitch_shift(3:end);
    std_targeted_pitch_shift(1) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1 plots averaged pitch shifts of targeted, same type and different
%type syllables combined across birds. Values are averages of mean pitch
%shifts of individual syllables.

len = length(mean_targeted_pitch_shift);

if birds_to_use(1) == 19
    figure (1)
    hold all
    legend_str = {};
    errorbar((1:len)-1,mean_targeted_pitch_shift,std_targeted_pitch_shift,'r','lineWidth',2)
    errorbar((1:len)-0.9,mean_same_pitch_shift,std_same_pitch_shift,'g','lineWidth',2)
    errorbar((1:len)-1.1,mean_diff_pitch_shift,std_diff_pitch_shift,'b','lineWidth',2)
    scatter((1:len)-1,mean_targeted_pitch_shift,30,'k','fill')
    scatter((1:len)-0.9,mean_same_pitch_shift,30,'k','fill')
    scatter((1:len)-1.1,mean_diff_pitch_shift,30,'k','fill')
    legend_str = [legend_str; 'Targeted Syllable';'Same type Syl';'Diff type Syl'];
    legend(legend_str)
    
    if plot_baseline
        xlim([0 (num_wn_days + num_base_days + 1)]) % make xlim one larger than number of elements (note that num_sets is one more than the # of files plotted for WN because the first file was just combined baseline) 
        set(gca,'Xtick',[1:(num_wn_days + num_base_days)]) % set xticks from 1:total-#-of-plotted-points
        xticklabel = {};
        for i=1:num_base_days
            xticklabel = [xticklabel; 'B ' num2str(i)];
        end
        for i=1:num_wn_days % start @ 2 because first file loaded in top_folder is just combined baseline days fiile, not a Wn day file
           xticklabel = [xticklabel; 'WN ' num2str(i)];
        end
        set(gca,'Xticklabel', xticklabel);
    else
        xlim([-0.2 len-0.8])
        ylim([-0.4 0.6])
    end

    % 0 cents line
     xl = xlim();
     plot(xlim, [0 0], 'k--');
    if plot_baseline
    % Baseline WN division
         yl = ylim();
         plot([3.5 3.5], yl, 'k--');
         ylabel('Adaptive Pitch Shift from Baseline in semitones')
    end

else
    figure (1)
    hold all
    legend_str = {};
    errorbar((1:len)-1,mean_targeted_pitch_shift,std_targeted_pitch_shift,'r','lineWidth',2)
    errorbar((1:len)-0.9,mean_same_pitch_shift,std_same_pitch_shift,'g','lineWidth',2)
    errorbar((1:len)-1.1,mean_diff_pitch_shift,std_diff_pitch_shift,'b','lineWidth',2)
    scatter((1:len)-1,mean_targeted_pitch_shift,30,'k','fill')
    scatter((1:len)-0.9,mean_same_pitch_shift,30,'k','fill')
    scatter((1:len)-1.1,mean_diff_pitch_shift,30,'k','fill')
    legend_str = [legend_str; 'Targeted Syllable';'Same type Syl';'Diff type Syl'];
    legend(legend_str)
    
    if plot_baseline
        xlim([0 (num_wn_days + num_base_days + 1)]) % make xlim one larger than number of elements (note that num_sets is one more than the # of files plotted for WN because the first file was just combined baseline) 
        set(gca,'Xtick',[1:(num_wn_days + num_base_days)]) % set xticks from 1:total-#-of-plotted-points
        xticklabel = {};
        for i=1:num_base_days
            xticklabel = [xticklabel; 'B ' num2str(i)];
        end
        for i=1:num_wn_days % start @ 2 because first file loaded in top_folder is just combined baseline days fiile, not a Wn day file
           xticklabel = [xticklabel; 'WN ' num2str(i)];
        end
        set(gca,'Xticklabel', xticklabel);
    else
        xlim([-0.2 len-0.8])
        ylim([-0.3 0.7])
    end

    % 0 cents line
     xl = xlim();
     plot(xlim, [0 0], 'k--');
    if plot_baseline
        % Baseline WN division
         yl = ylim();
         plot([3.5 3.5], yl, 'k--');
         ylabel('Adaptive Pitch Shift from Baseline in semitones')
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following section of code is used to plot the above plot as it would
%have been if calculated simply by using Miterations as in previous papers.
%Use this section by toggling the plot_non_bootstrapped variable set at the
%beginning of the code.

if plot_non_bootstrapped
    mean_targeted_pitch_shift_nb = zeros(num_base_days + num_wn_days,1);
    mean_same_pitch_shift_nb = zeros(num_base_days + num_wn_days,1);
    mean_diff_pitch_shift_nb = zeros(num_base_days + num_wn_days,1);
    std_targeted_pitch_shift_nb = zeros(num_base_days + num_wn_days,1);
    std_same_pitch_shift_nb = zeros(num_base_days + num_wn_days,1);
    std_diff_pitch_shift_nb = zeros(num_base_days + num_wn_days,1);
    
    for i = 1:num_days
        target_temp = [];
        same_temp = [];
        diff_temp = [];
        
        for j = 1:length(birds_to_use)
            for ind = 1:length(P(j).target_syl)
                temp = P(j).data{i,P(j).target_syl(ind)};
                if size(temp,1) == 1
                    temp = temp';
                end
                target_temp = [target_temp; temp];
            end
            for ind = 1:length(P(j).same_syls)
                temp = P(j).data{i,P(j).same_syls(ind)};
                if size(temp,1) == 1
                    temp = temp';
                end
                same_temp = [same_temp; temp];
            end
            for ind = 1:length(P(j).diff_syls)
                temp = P(j).data{i,P(j).diff_syls(ind)};
                if size(temp,1) == 1
                    temp = temp';
                end
                diff_temp = [diff_temp; temp];
            end
        end
        
        mean_targeted_pitch_shift_nb(i) = mean(target_temp);
        mean_same_pitch_shift_nb(i) = mean(same_temp);
        mean_diff_pitch_shift_nb(i) = mean(diff_temp);
        std_targeted_pitch_shift_nb(i) = std(target_temp)/length(target_temp);
        std_same_pitch_shift_nb(i) = std(same_temp)/length(same_temp);
        std_diff_pitch_shift_nb(i) = std(diff_temp)/length(diff_temp);
    end
    
    if plot_baseline == 0
        mean_same_pitch_shift_nb = mean_same_pitch_shift_nb(3:end);
        mean_same_pitch_shift_nb(1) = 0;
        std_same_pitch_shift_nb = std_same_pitch_shift_nb(3:end);
        std_same_pitch_shift_nb(1) = 0;
        mean_diff_pitch_shift_nb = mean_diff_pitch_shift_nb(3:end);
        mean_diff_pitch_shift_nb(1) = 0;
        std_diff_pitch_shift_nb = std_diff_pitch_shift_nb(3:end);
        std_diff_pitch_shift_nb(1) = 0;
        mean_targeted_pitch_shift_nb = mean_targeted_pitch_shift_nb(3:end);
        mean_targeted_pitch_shift_nb(1) = 0;
        std_targeted_pitch_shift_nb = std_targeted_pitch_shift_nb(3:end);
        std_targeted_pitch_shift_nb(1) = 0;
    end
    
    figure (2)
    hold all
    errorbar((1:len)-1,mean_targeted_pitch_shift_nb,std_targeted_pitch_shift_nb,'r','lineWidth',2,'DisplayName','Targeted Syllable')
    errorbar((1:len)-0.9,mean_same_pitch_shift_nb,std_same_pitch_shift_nb,'g','lineWidth',2,'DisplayName','Same type Syllable')
    errorbar((1:len)-1.1,mean_diff_pitch_shift_nb,std_diff_pitch_shift_nb,'b','lineWidth',2,'DisplayName','Diff type Syllable')
%     scatter(1:len,mean_targeted_pitch_shift_nb,30,'k','fill')
%     scatter((1:len)+0.1,mean_same_pitch_shift_nb,30,'k','fill')
%     scatter((1:len)-0.1,mean_diff_pitch_shift_nb,30,'k','fill')
    legend()

    if plot_baseline
        xlim([0 (num_wn_days + num_base_days + 1)]) % make xlim one larger than number of elements (note that num_sets is one more than the # of files plotted for WN because the first file was just combined baseline) 
        set(gca,'Xtick',[1:(num_wn_days + num_base_days)]) % set xticks from 1:total-#-of-plotted-points
        xticklabel = {};
        for i=1:num_base_days
            xticklabel = [xticklabel; 'B ' num2str(i)];
        end
        for i=1:num_wn_days % start @ 2 because first file loaded in top_folder is just combined baseline days fiile, not a Wn day file
           xticklabel = [xticklabel; 'WN ' num2str(i)];
        end
        set(gca,'Xticklabel', xticklabel);
    else
        xlim([-0.2 len-0.8])
    end

    % 0 cents line
     xl = xlim();
     plot(xlim, [0 0], 'k--');
    if plot_baseline
        % Baseline WN division
         yl = ylim();
         plot([3.5 3.5], yl, 'k--');
         ylabel('Adaptive Pitch Shift from Baseline in semitones')
    end
end

%%
%We additionally need a section to compute the bootstrap without the
%hierarchical structure particularly for the WN data. That will be done
%here. We will also compute statistics from this subsection:

rng(1); %For reproducibility of stats.

std_targeted_pitch_shift_wrong_boot = zeros(num_base_days + num_wn_days,1);
std_same_pitch_shift_wrong_boot = zeros(num_base_days + num_wn_days,1);
std_diff_pitch_shift_wrong_boot = zeros(num_base_days + num_wn_days,1);

nboot = 1000;

for i = 1:num_days
    temp2 = [];
    for j = 1:length(birds_to_use)
        if birds_to_use(1) == 19
            temp = horzcat(P(j).data{i,P(j).target_syl});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(j).data{i,P(j).target_syl});
            temp2 = vertcat(temp2,temp);
        end
    end
    bootstats1 = bootstrp(nboot,@mean,temp2);
    std_targeted_pitch_shift_wrong_boot(i) = std(bootstats1);

    temp2 = [];
    for j = 1:length(same_birds)
        if birds_to_use(1) == 19
            temp = horzcat(P(same_birds(j)).data{i,P(same_birds(j)).same_syls});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(same_birds(j)).data{i,P(same_birds(j)).same_syls});
            temp2 = vertcat(temp2,temp);
        end
    end
    bootstats2 = bootstrp(nboot,@mean,temp2);
    std_same_pitch_shift_wrong_boot(i) = std(bootstats2);

    temp2 = [];
    for j = 1:length(diff_birds)
        if birds_to_use(1) == 19
            temp = horzcat(P(diff_birds(j)).data{i,P(diff_birds(j)).diff_syls});
            temp2 = horzcat(temp2,temp);
        else
            temp = vertcat(P(diff_birds(j)).data{i,P(diff_birds(j)).diff_syls});
            temp2 = vertcat(temp2,temp);
        end
    end
    bootstats3 = bootstrp(nboot,@mean,temp2);
    std_diff_pitch_shift_wrong_boot(i) = std(bootstats3);

end

if plot_baseline == 0
    std_same_pitch_shift_wrong_boot = std_same_pitch_shift_wrong_boot(3:end);
    std_same_pitch_shift_wrong_boot(1) = 0;
    std_diff_pitch_shift_wrong_boot = std_diff_pitch_shift_wrong_boot(3:end);
    std_diff_pitch_shift_wrong_boot(1) = 0;
    std_targeted_pitch_shift_wrong_boot = std_targeted_pitch_shift_wrong_boot(3:end);
    std_targeted_pitch_shift_wrong_boot(1) = 0;
end

len = length(mean_targeted_pitch_shift);

if len ~= length(std_same_pitch_shift_wrong_boot)
    mean_same_pitch_shift = mean_same_pitch_shift(3:end);
    mean_same_pitch_shift(1) = 0;
    mean_diff_pitch_shift = mean_diff_pitch_shift(3:end);
    mean_diff_pitch_shift(1) = 0;
    mean_targeted_pitch_shift = mean_targeted_pitch_shift(3:end);
    mean_targeted_pitch_shift(1) = 0;
end

len = length(mean_targeted_pitch_shift);

figure (3)
hold all
legend_str = {};
errorbar((1:len)-1,mean_targeted_pitch_shift,std_targeted_pitch_shift_wrong_boot,'r','lineWidth',2)
errorbar((1:len)-0.9,mean_same_pitch_shift,std_same_pitch_shift_wrong_boot,'g','lineWidth',2)
errorbar((1:len)-1.1,mean_diff_pitch_shift,std_diff_pitch_shift_wrong_boot,'b','lineWidth',2)
scatter((1:len)-1,mean_targeted_pitch_shift,30,'k','fill')
scatter((1:len)-0.9,mean_same_pitch_shift,30,'k','fill')
scatter((1:len)-1.1,mean_diff_pitch_shift,30,'k','fill')
legend_str = [legend_str; 'Targeted Syllable';'Same type Syl';'Diff type Syl'];
legend(legend_str)

if plot_baseline
    xlim([0 (num_wn_days + num_base_days + 1)]) % make xlim one larger than number of elements (note that num_sets is one more than the # of files plotted for WN because the first file was just combined baseline) 
    set(gca,'Xtick',[1:(num_wn_days + num_base_days)]) % set xticks from 1:total-#-of-plotted-points
    xticklabel = {};
    for i=1:num_base_days
        xticklabel = [xticklabel; 'B ' num2str(i)];
    end
    for i=1:num_wn_days % start @ 2 because first file loaded in top_folder is just combined baseline days fiile, not a Wn day file
       xticklabel = [xticklabel; 'WN ' num2str(i)];
    end
    set(gca,'Xticklabel', xticklabel);
else
    xlim([-0.2 len-0.8])
    ylim([-0.3 0.7])
end

% 0 cents line
 xl = xlim();
 plot(xlim, [0 0], 'k--');
if plot_baseline
% Baseline WN division
     yl = ylim();
     plot([3.5 3.5], yl, 'k--');
     ylabel('Adaptive Pitch Shift from Baseline in semitones')
end

sum(bootstats1>=0)/length(bootstats1)
sum(bootstats2>=0)/length(bootstats2)
sum(bootstats3>=0)/length(bootstats3)


%%
%We should also produce a figure displaying the means and SEMs for both
%cases when we use the summarized method. That will be done in this
%section.

targeted_syl_means = zeros(num_wn_days+1,1);
same_syl_means = zeros(num_wn_days+1,1);
diff_syl_means = zeros(num_wn_days+1,1);

targeted_syl_sem = zeros(num_wn_days+1,1);
same_syl_sem = zeros(num_wn_days+1,1);
diff_syl_sem = zeros(num_wn_days+1,1);

for i = (num_base_days+1):num_days
    temp2 = [];
    for j = 1:length(birds_to_use)
        for syl = P(j).target_syl
            
            temp = mean(P(j).data{i,syl});
            temp2 = [temp2; temp];
                       
        end
    end
    targeted_syl_means(i-num_base_days+1) = mean(temp2);
    targeted_syl_sem(i-num_base_days+1) = std(temp2)/sqrt(length(temp2));
    
    temp2 = [];
    for j = 1:length(same_birds)
        for syl = P(same_birds(j)).same_syls
            
            temp = mean(P(same_birds(j)).data{i,syl});
            temp2 = [temp2; temp];
                       
        end
    end
    same_syl_means(i-num_base_days+1) = mean(temp2);
    same_syl_sem(i-num_base_days+1) = std(temp2)/sqrt(length(temp2));
    
    temp2 = [];
    for j = 1:length(diff_birds)
        for syl = P(diff_birds(j)).diff_syls
            
            temp = mean(P(diff_birds(j)).data{i,syl});
            temp2 = [temp2; temp];
                       
        end
    end
    diff_syl_means(i-num_base_days+1) = mean(temp2);
    diff_syl_sem(i-num_base_days+1) = std(temp2)/sqrt(length(temp2));

end

len = length(targeted_syl_means);
legend_str = {};
figure()
hold all
errorbar((1:len)-1,targeted_syl_means,targeted_syl_sem,'r','lineWidth',2)
errorbar((1:len)-0.9,same_syl_means,same_syl_sem,'g','lineWidth',2)
errorbar((1:len)-1.1,diff_syl_means,diff_syl_sem,'b','lineWidth',2)
scatter((1:len)-1,targeted_syl_means,30,'k','fill')
scatter((1:len)-0.9,same_syl_means,30,'k','fill')
scatter((1:len)-1.1,diff_syl_means,30,'k','fill')
legend_str = [legend_str; 'Targeted Syllable';'Same type Syl';'Diff type Syl'];
legend(legend_str)

if plot_baseline
    xlim([0 (num_wn_days + num_base_days + 1)]) % make xlim one larger than number of elements (note that num_sets is one more than the # of files plotted for WN because the first file was just combined baseline) 
    set(gca,'Xtick',[1:(num_wn_days + num_base_days)]) % set xticks from 1:total-#-of-plotted-points
    xticklabel = {};
    for i=1:num_base_days
        xticklabel = [xticklabel; 'B ' num2str(i)];
    end
    for i=1:num_wn_days % start @ 2 because first file loaded in top_folder is just combined baseline days fiile, not a Wn day file
       xticklabel = [xticklabel; 'WN ' num2str(i)];
    end
    set(gca,'Xticklabel', xticklabel);
else
    xlim([-0.2 len-0.8])
    ylim([-0.3 0.7])
end

% 0 cents line
 xl = xlim();
 plot(xlim, [0 0], 'k--');
if plot_baseline
% Baseline WN division
     yl = ylim();
     plot([3.5 3.5], yl, 'k--');
     ylabel('Adaptive Pitch Shift from Baseline in semitones')
end



%%
%We also need a section to compute stats using the Traditional and
%Summarized metrics. We will have to report them.

%First, we will do the Traditional method. For this, we need only condense
%all the data for the relevant group and perform a t-test comparing for
%difference to zero.

all_targeted_syllables = [];
all_same_syllables = [];
all_diff_syllables = [];

if wn
    subtract_days = 1;
else
    subtract_days = 2;
end

for i = (num_days-subtract_days):num_days
    for j = 1:length(birds_to_use)
        if birds_to_use(1) == 19
            temp = horzcat(P(j).data{i,P(j).target_syl});
            all_targeted_syllables = horzcat(all_targeted_syllables,temp);
        else
            temp = vertcat(P(j).data{i,P(j).target_syl});
            all_targeted_syllables = vertcat(all_targeted_syllables,temp);
        end
    end
    
    for j = 1:length(same_birds)
        if birds_to_use(1) == 19
            temp = horzcat(P(same_birds(j)).data{i,P(same_birds(j)).same_syls});
            all_same_syllables = horzcat(all_same_syllables,temp);
        else
            temp = vertcat(P(same_birds(j)).data{i,P(same_birds(j)).same_syls});
            all_same_syllables = vertcat(all_same_syllables,temp);
        end
    end
    
    for j = 1:length(diff_birds)
        if birds_to_use(1) == 19
            temp = horzcat(P(diff_birds(j)).data{i,P(diff_birds(j)).diff_syls});
            all_diff_syllables = horzcat(all_diff_syllables,temp);
        else
            temp = vertcat(P(diff_birds(j)).data{i,P(diff_birds(j)).diff_syls});
            all_diff_syllables = vertcat(all_diff_syllables,temp);
        end
    end
    
end

summary_stats_trad = zeros(1,6);
summary_stats_trad(1) = mean(all_targeted_syllables);
summary_stats_trad(2) = std(all_targeted_syllables)/sqrt(length(all_targeted_syllables));
summary_stats_trad(3) = mean(all_same_syllables);
summary_stats_trad(4) = std(all_same_syllables)/sqrt(length(all_same_syllables));
summary_stats_trad(5) = mean(all_diff_syllables);
summary_stats_trad(6) = std(all_diff_syllables)/sqrt(length(all_diff_syllables));

[~,p1_trad,~,stats1_trad] = ttest(all_targeted_syllables)
[~,p2_trad,~,stats2_trad] = ttest(all_same_syllables)
[~,p3_trad,~,stats3_trad] = ttest(all_diff_syllables)

%Now for the summarized statistics. Here we need to combine all the
%syllables for a given type and extract only the mean. All the means are
%then saved and used for analysis.

targeted_syllable_means = [];
same_syllable_means = [];
diff_syllable_means = [];

for j = 1:length(birds_to_use)
    for syl = P(j).target_syl
        if birds_to_use(1) == 19
            temp = horzcat(P(j).data{end-2:end,syl});
            targeted_syllable_means = [targeted_syllable_means; mean(temp)];
        else
            temp = vertcat(P(j).data{end-1:end,syl});
            targeted_syllable_means = [targeted_syllable_means; mean(temp)];
        end
    end
end

for j = 1:length(same_birds)
    for syl = P(same_birds(j)).same_syls
        if birds_to_use(1) == 19
            temp = horzcat(P(same_birds(j)).data{end-2:end,syl});
            same_syllable_means = [same_syllable_means; mean(temp)];
        else
            temp = vertcat(P(same_birds(j)).data{end-1:end,syl});
            same_syllable_means = [same_syllable_means; mean(temp)];
        end
    end
end

for j = 1:length(diff_birds)
    for syl = P(diff_birds(j)).diff_syls
        if birds_to_use(1) == 19
            temp = horzcat(P(diff_birds(j)).data{end-2:end,syl});
            diff_syllable_means = [diff_syllable_means; mean(temp)];
        else
            temp = vertcat(P(diff_birds(j)).data{end-1:end,syl});
            diff_syllable_means = [diff_syllable_means; mean(temp)];
        end
    end
end

summary_stats_summ = zeros(1,6);
summary_stats_summ(1) = mean(targeted_syllable_means);
summary_stats_summ(2) = std(targeted_syllable_means)/sqrt(length(targeted_syllable_means));
summary_stats_summ(3) = mean(same_syllable_means);
summary_stats_summ(4) = std(same_syllable_means)/sqrt(length(same_syllable_means));
summary_stats_summ(5) = mean(diff_syllable_means);
summary_stats_summ(6) = std(diff_syllable_means)/sqrt(length(diff_syllable_means));

[~,p1_summ,~,stats1_summ] = ttest(targeted_syllable_means)
[~,p2_summ,~,stats2_summ] = ttest(same_syllable_means)
[~,p3_summ,~,stats3_summ] = ttest(diff_syllable_means)

%%
%This section is used to compute stats. We set the random number generator
%in this section so as to ensure reproducibility.

rng(1);

%Bootstrapping part
nboot = 300; %No of times to resample for bootstrapping.
max_syls = 16;

bootstrapping_matrix = zeros(nboot,(num_days+1),max_syls,numel(birds_to_use));

for i = 1:numel(birds_to_use)   %Over all birds
    for j = 1:size(P(i).data,2) %Over number of syllables for each bird
        for k = 1:num_days      %Over all days
            temp = P(i).data_hz{k,j}; %Remember to take actual Hz values for this step.
            if isempty(temp)
                bootstrapping_matrix(:,k,j,i) = NaN;
                continue;
            end
            bootstrapping_matrix(:,k,j,i) = datasample(temp,nboot);
        end
        for n = 1:nboot
            if birds_to_use(1) == 19 %Average over last 3 days for headphones
                bootstrapping_matrix(n,end,j,i) = nanmean(bootstrapping_matrix(n,(end-3):(end-1),j,i));
            else %Average over last two days for WN
                bootstrapping_matrix(n,end,j,i) = nanmean(bootstrapping_matrix(n,(end-2):(end-1),j,i));
            end
            temp_mean = nanmean(bootstrapping_matrix(n,1:3,j,i)); %Normalize each row to the mean of the 3 baseline days of that row.
            if birds_to_use(1) == 19
                bootstrapping_matrix(n,:,j,i) = -1*P(i).shift_direction*12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            else
                bootstrapping_matrix(n,:,j,i) = P(i).shift_direction*12*log2(bootstrapping_matrix(n,:,j,i)/temp_mean);
            end
        end
    end
end

nboot2 = 10000; %Repeating 10000 times since this is the distribution used to calculate stats.
bootstats1 = zeros(nboot2,1);
bootstats2 = zeros(nboot2,1);
bootstats3 = zeros(nboot2,1);
for n = 1:nboot2
    
    temp_targ_birds = datasample(1:numel(birds_to_use),numel(birds_to_use));
    temp_data = [];
    for t = 1:length(temp_targ_birds)
        temp_syls = datasample(P(temp_targ_birds(t)).target_syl,length(P(temp_targ_birds(t)).target_syl));
        for s = 1:length(temp_syls)
            if birds_to_use(1) == 19
                num_iterations = length(horzcat(P(temp_targ_birds(t)).data_hz{(end-2):end,temp_syls(s)}));
            else
                num_iterations = length(vertcat(P(temp_targ_birds(t)).data_hz{(end-1):end,temp_syls(s)}));
            end
            if num_iterations == 0
                num_iterations = 1;
            end
            temp_pulls = datasample(1:nboot,num_iterations);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_targ_birds(t)));
        end
    end
    bootstats1(n,1) = nanmean(temp_data);
    
    temp_diff_birds = datasample(diff_birds,length(diff_birds));
    temp_data = [];
    for t = 1:length(temp_diff_birds)
        temp_syls = datasample(P(temp_diff_birds(t)).diff_syls,length(P(temp_diff_birds(t)).diff_syls));
        for s = 1:length(temp_syls)
            if birds_to_use(1) == 19
                num_iterations = length(horzcat(P(temp_diff_birds(t)).data_hz{(end-2):end,temp_syls(s)}));
            else
                num_iterations = length(vertcat(P(temp_diff_birds(t)).data_hz{(end-1):end,temp_syls(s)}));
            end
            if num_iterations == 0
                num_iterations = 1;
            end
            temp_pulls = datasample(1:nboot,num_iterations);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_diff_birds(t)));
        end
    end
    bootstats2(n,1) = nanmean(temp_data);
    
    temp_same_birds = datasample(same_birds,length(same_birds));
    temp_data = [];
    for t = 1:length(temp_same_birds)
        temp_syls = datasample(P(temp_same_birds(t)).same_syls,length(P(temp_same_birds(t)).same_syls));
        for s = 1:length(temp_syls)
            if birds_to_use(1) == 19
                num_iterations = length(horzcat(P(temp_same_birds(t)).data_hz{(end-2):end,temp_syls(s)}));
            else
                num_iterations = length(vertcat(P(temp_same_birds(t)).data_hz{(end-1):end,temp_syls(s)}));
            end
            if num_iterations == 0
                num_iterations = 1;
            end
            temp_pulls = datasample(1:nboot,num_iterations);
            temp_data = vertcat(temp_data,bootstrapping_matrix(temp_pulls,end,temp_syls(s),temp_same_birds(t)));
        end
    end
    bootstats3(n,1) = nanmean(temp_data);
end

sum(bootstats1>=0)/length(bootstats1)
sum(bootstats2>=0)/length(bootstats2)
sum(bootstats3>=0)/length(bootstats3)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We now use a similar bootstrapping strategy as defined above to calculate
%the mean and error bars for final syllable shifts based on position from
%the target syllable. This could get a little tricky due to the rarity of
%each type of syllable at different positions but we do have all the
%required information in the structure variable P.

nboot = 1000; %No of times to resample for bootstrapping.
lower_dist_limit = min(horzcat(P.syls_away_from_target));
upper_dist_limit = max(horzcat(P.syls_away_from_target));
distance_graphs_x = lower_dist_limit:upper_dist_limit;
distance_graphs_same_pitch = zeros(size(distance_graphs_x));
distance_graphs_diff_pitch = zeros(size(distance_graphs_x));
distance_graphs_same_std = zeros(size(distance_graphs_x));
distance_graphs_diff_std = zeros(size(distance_graphs_x));

for i = 1:length(distance_graphs_x)
    
    same_birds = [];
    diff_birds = [];
    for j = 1:num_birds
        if sum(find(P(j).syls_away_from_target == distance_graphs_x(i)))
            ind = find(P(j).syls_away_from_target == distance_graphs_x(i));
            if sum(find(P(j).same_syls == ind(1)))
                for dei = 1:length(ind)
                    same_birds = [same_birds; [j ind(dei)]];
                end
            elseif sum(find(P(j).diff_syls == ind(1)))
                for dei = 1:length(ind)
                    diff_birds = [diff_birds; [j ind(dei)]];
                end
            end
        end
    end
    
    bootstats_same = zeros(nboot,1);
    if ~(isempty(same_birds))
        for n = 1:nboot

            temp_birds = datasample(1:size(same_birds,1),size(same_birds,1));
            temp_data = [];
            for t = 1:length(temp_birds)
                temp = P(same_birds(temp_birds(t),1)).data{end,same_birds(temp_birds(t),2)};
                temp2 = datasample(temp,length(temp));
                temp_data = [temp_data; temp2];
            end
            bootstats_same(n) = mean(temp_data);

        end
        distance_graphs_same_pitch(i) = mean(bootstats_same);
        distance_graphs_same_std(i) = 1.96*std(bootstats_same);
    else
        distance_graphs_same_pitch(i) = NaN;
        distance_graphs_same_std(i) = NaN;
    end
    
    bootstats_diff = zeros(nboot,1);
    if ~(isempty(diff_birds))
        for n = 1:nboot

            temp_birds = datasample(1:size(diff_birds,1),size(diff_birds,1));
            temp_data = [];
            for t = 1:length(temp_birds)
                temp = P(diff_birds(temp_birds(t),1)).data{end,diff_birds(temp_birds(t),2)};
                temp2 = datasample(temp,length(temp));
                temp_data = [temp_data; temp2];
            end
            bootstats_diff(n) = mean(temp_data);

        end
        distance_graphs_diff_pitch(i) = mean(bootstats_diff);
        distance_graphs_diff_std(i) = 1.96*std(bootstats_diff);
    else
        distance_graphs_diff_pitch(i) = NaN;
        distance_graphs_diff_std(i) = NaN;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figures 2 and 3 plot how same type and different type syllables vary as a
%function of distance from the target syllable across birds.

figure (2)
errorbar(distance_graphs_x,distance_graphs_same_pitch,distance_graphs_same_std,'g','LineWidth',2)
hold on
scatter(distance_graphs_x,distance_graphs_same_pitch,30,'k','fill')
errorbar(0,mean_targeted_pitch_shift(end),std_targeted_pitch_shift(end),'r','lineWidth',2)
scatter(0,mean_targeted_pitch_shift(end),10,'k','fill')
plot([lower_dist_limit upper_dist_limit],[0 0], 'k', 'lineWidth', 2)
xlabel('Sequential Distance from target syllable')
ylabel('Adaptive pitch change from baseline in semitones')
title('Same type syllables with syllable distance')
% print(f5,'-depsc',['fig5.eps']);

figure (3)
errorbar(distance_graphs_x,distance_graphs_diff_pitch,distance_graphs_diff_std,'b','LineWidth',2)
hold on
scatter(distance_graphs_x,distance_graphs_diff_pitch,30,'k','fill')
errorbar(0,mean_targeted_pitch_shift(end),std_targeted_pitch_shift(end),'r','lineWidth',2)
scatter(0,mean_targeted_pitch_shift(end),10,'k','fill')
plot([lower_dist_limit upper_dist_limit],[0 0], 'k', 'lineWidth', 2)
xlabel('Sequential Distance from target syllable')
ylabel('Adaptive pitch change from baseline in semitones')
title('Different type syllables with syllable distance')
% print(f6,'-depsc',['fig6.eps']);

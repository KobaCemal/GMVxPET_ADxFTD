clear

%% Get the necessary tables
subjects_fs_data=readtable('/home/koba/Desktop/kramer/adni/itsavibeNature2/UCSFFSX7_05May2025.csv');
region_codes=readtable('/home/koba/Desktop/kramer/adni/itsavibeNature2/full_region_table_with_codes.csv');
demographics=readtable('/home/koba/Desktop/kramer/adni/itsavibeNature2/ADNIMERGE/ADNIMERGE_05May2025 (4).csv');
petdata_full=readtable('/home/koba/Desktop/kramer/JuSpace-JuSpace_v1.5/JuSpace_v1.5/txtfiles/parcellated_petfiles.csv');
pet_of_interest=[1 4 6 10 11 13 14 15 18 20 24 25 30];
petdata=petdata_full(pet_of_interest,:);
petdata_mat=table2array(petdata(:,2:end));

% add Combat toolbox to the path
addpath(genpath('/home/koba/Desktop/kramer/adni/itsavibeNature2/ComBatHarmonization-master/'))
addpath(genpath('/home/koba/Desktop/kramer/adni/itsavibeNature2/ENIGMA'))

%% Prepare the Freesurfer data as we had in the previous study (83 regions)
% Extract the 'New Code' column
new_codes = region_codes.newcode;

% Remove entries that are "NOT FOUND" or empty
valid_codes = new_codes(~strcmp(new_codes, 'NOT FOUND') & ~cellfun(@isempty, new_codes));

% Keep columns in the same order as 'valid_codes'
[~, ia, ib] = intersect(valid_codes, subjects_fs_data.Properties.VariableNames, 'stable');

% Subset the data table
subjects_fs_data_filtered = subjects_fs_data(:, valid_codes(ia));
subjects_fs_data_filtered.BrainStem=subjects_fs_data_filtered.ST76SV + subjects_fs_data_filtered.ST17SV; % Brain stem is given in two hemispheres. Summing them up to be consistent with the previous analyses. 
subjects_fs_data_filtered.ST76SV=[];
subjects_fs_data_filtered.ST17SV=[];

% Number of Nans and replacement with mean 
sum(sum(isnan(table2array(subjects_fs_data_filtered)))) % 2044  -- 0.4 percent of the total values
max(sum(isnan(table2array(subjects_fs_data_filtered)))) % max number of nans -- 65  
subjects_fs_data_filtered_mat=table2array(subjects_fs_data_filtered); % get an array version for ease of use

% Replace missing values 
for i = 1:width(subjects_fs_data_filtered_mat)
    column = subjects_fs_data_filtered_mat(:, i);

    if isnumeric(column)
        % Compute column mean ignoring NaNs
        col_mean = mean(column, 'omitnan');

        % Replace NaNs with the column mean
        column(isnan(column)) = col_mean;

        % Assign the cleaned column back to the table
        subjects_fs_data_filtered_mat(:, i)= column;
    end
end

cortical_index = strcmp(region_codes.structure,'cortex'); % must get the cortical and subcortical region index from the region list for visualization
subcortical_index=[find(strcmp(region_codes.label,'accumbensarea') & strcmp(region_codes.hemisphere, 'L'))...
    find(strcmp(region_codes.label,'amygdala') & strcmp(region_codes.hemisphere, 'L'))...
    find(strcmp(region_codes.label,'caudate') & strcmp(region_codes.hemisphere, 'L'))...
    find(strcmp(region_codes.label,'hippocampus') & strcmp(region_codes.hemisphere, 'L'))...
    find(strcmp(region_codes.label,'pallidum') & strcmp(region_codes.hemisphere, 'L'))...
    find(strcmp(region_codes.label,'putamen') & strcmp(region_codes.hemisphere, 'L'))...
    find(strcmp(region_codes.label,'thalamusproper') & strcmp(region_codes.hemisphere, 'L'))...
    find(strcmp(region_codes.label,'accumbensarea') & strcmp(region_codes.hemisphere, 'R'))...
    find(strcmp(region_codes.label,'amygdala') & strcmp(region_codes.hemisphere, 'R'))...
    find(strcmp(region_codes.label,'caudate') & strcmp(region_codes.hemisphere, 'R'))...
    find(strcmp(region_codes.label,'hippocampus') & strcmp(region_codes.hemisphere, 'R'))...
    find(strcmp(region_codes.label,'pallidum') & strcmp(region_codes.hemisphere, 'R'))...
    find(strcmp(region_codes.label,'putamen') & strcmp(region_codes.hemisphere, 'R'))...
    find(strcmp(region_codes.label,'thalamusproper') & strcmp(region_codes.hemisphere, 'R'))]';
%% Get the demographics/metadata during the first session of the subject

% Create an empty row for replacing the subjects with missing data
empty_row = demographics(1, :);  
for var = 1:width(empty_row)
    val = empty_row{1, var};

    if isnumeric(val) || islogical(val)
        empty_row{1, var} = NaN;
    elseif isstring(val)
        empty_row{1, var} = missing;
    elseif iscellstr(val)
        empty_row{1, var} = {''};
    elseif iscategorical(val)
        empty_row{1, var} = categorical(missing);
    else
        empty_row{1, var} = missing;  % fallback
    end
end

% Every subjct has multiple sessions. Choose the one with 'bl' code in
% VISCODE column
subject_demographics_single_session=repmat(empty_row,size(subjects_fs_data_filtered,1),1);
missingsubjects=nan(size(subjects_fs_data_filtered,1),1);

for i=1:size(subjects_fs_data_filtered,1)
    if sum(contains(demographics.PTID,subjects_fs_data.PTID{i}))>0 %if the data has its metadata
        demo=demographics(contains(demographics.PTID,subjects_fs_data.PTID{i}),:); % get the metadata across sessions
        if sum(strcmp(demo.VISCODE,'bl')) == 1
            subject_demographics_single_session(i,:)=demo(strcmp(demo.VISCODE,'bl'),:); % choose the 'bl' session. if not, mark it as missing
            missingsubjects(i)=0;
        else
            %subject_demographics_single_session(i,:)=empty_row;
            missingsubjects(i)=1;
        end
    else
        missingsubjects(i)=1;
    end
end

% Remove the subjects that dont have any metadata from both pools
sum(missingsubjects) % 339 subjects are removed here

subject_demographics_single_session(logical(missingsubjects'),:)=[];
subjects_fs_data_filtered_mat(logical(missingsubjects'),:)=[];

missing_demographics= ismissing(subject_demographics_single_session.AGE)...
    + ismissing(subject_demographics_single_session.PTGENDER)...
    + ismissing(subject_demographics_single_session.ICV)...
    + ismissing(subject_demographics_single_session.SITE) ...
    + ismissing(subject_demographics_single_session.DX_bl);
sum(missing_demographics) % 168 subjects with at least one missing demographics

subject_demographics_single_session(missing_demographics>0,:)=[];
subjects_fs_data_filtered_mat(missing_demographics>0,:)=[];

% get age, sex, and intracranial volume (ICV -- confound regressors we used in the previous study)
% -- in addition, we need to get the site info, as there are dat afrom
% multiple scanners
age=subject_demographics_single_session.AGE;
sex=subject_demographics_single_session.PTGENDER;
sex_dummy=grp2idx(sex)-1;
icv=subject_demographics_single_session.ICV;
site=subject_demographics_single_session.SITE;
diagnosis=subject_demographics_single_session.DX_bl; % there is another DX column with different labeling

labels=unique(subject_demographics_single_session.DX_bl); % unique labels: AD, CN, EMCI, LMCI, SMC

%% Ensure each site has at least 10 observations
site=subject_demographics_single_session.SITE;
diagnosis=subject_demographics_single_session.DX_bl; % there is another DX column with different labeling

% Combine into a table for processing
T = table(site, diagnosis);

% Step 1: Get unique site-diagnosis combinations and their counts
[G, site_list, dx_list] = findgroups(T.site, T.diagnosis);
counts = splitapply(@numel, T.site, G);

% Step 2: Build a table of counts
count_table = table(site_list, dx_list, counts);

% Step 3: Find sites that have at least 5 for *each* diagnosis
min_required = 5;
unique_sites = unique(site);

valid_sites = [];
for i = 1:length(unique_sites)
    this_site = unique_sites(i);
    rows = count_table.site_list == this_site;

    % Get all diagnoses for this site
    site_dx = count_table(rows, :);

    % Check if all 5 categories are present with enough subjects
    dx_names = {'AD', 'CN', 'EMCI', 'LMCI', 'SMC'};
    ok = true;

    for j = 1:length(dx_names)
        match = strcmp(site_dx.dx_list, dx_names{j});
        if sum(match) == 0 || site_dx.counts(match) < min_required
            ok = false;
            break;
        end
    end

    if ok
        valid_sites(end+1) = this_site;
    end
end

% Step 4: Filter original table to include only valid sites
keep_idx = ismember(site, valid_sites);
T_reduced = T(keep_idx, :);
subjects_data_reduced=subjects_fs_data_filtered_mat(keep_idx,:);
demographics_reduced=subject_demographics_single_session(keep_idx,:);
%% Initial analyses on demographics 
age=demographics_reduced.AGE;
sex=demographics_reduced.PTGENDER;
sex_dummy=grp2idx(sex)-1;
icv=demographics_reduced.ICV;
site=demographics_reduced.SITE;
diagnosis=demographics_reduced.DX_bl; % there is another DX column with different labeling
tau=demographics_reduced.TAU;
ptau=demographics_reduced.PTAU;
mmse=demographics_reduced.MMSE;
MOCA=demographics_reduced.MOCA;
adascog=demographics_reduced.ADAS11;
AV45_bl=demographics_reduced.AV45_bl;


labels=unique(demographics_reduced.DX_bl);
sum(strcmp(sex,'Male'))
sum(strcmp(sex,'Female'))
mean(age)
std(age)

[TABLE,CHI2,P,LABELS] = crosstab(sex,diagnosis) % sex

[P,ANOVATAB,STATS] =anova1(age,diagnosis); % age 
multcompare(STATS) % The groups have different sex  distribution 
mean(age(strcmp(diagnosis,'AD')))
std(age(strcmp(diagnosis,'AD'))) % repeat with other groups for methods section

[P,ANOVATAB,STATS] =anova1(icv,diagnosis); % icv 
multcompare(STATS) % The groups have different icv distribution 
mean(icv(strcmp(diagnosis,'AD')))
std(icv(strcmp(diagnosis,'AD')))

[P,ANOVATAB,STATS] =anova1(adascog(~isnan(adascog)),diagnosis(~isnan(adascog))); % icv 
multcompare(STATS) % The groups have different icv distribution 
mean(icv(strcmp(diagnosis,'AD')))
std(icv(strcmp(diagnosis,'AD')))

%% Remove the nuisance variables from GM data -- using combat for site effects

gmv_data=subjects_data_reduced;
batch=site';
confounds=[age sex_dummy icv];
diagnosis_nonan=diagnosis;

residuals=combat(gmv_data', batch, confounds, 1);


% Before/after combat
% Compute the mean GMV per subject (row-wise)
gmv_mean = mean(gmv_data, 2, 'omitnan');  % mean across columns for each row
boxplot(gmv_mean, batch)
xlabel('Site')
ylabel('Mean GMV')
ylim([5000 9000])
title('Mean GMV by Site / Before Combat')
figure
gmv_mean = mean(residuals, 'omitnan');  % mean across columns for each row
boxplot(gmv_mean, batch)
xlabel('Site')
ylabel('Mean GMV')
ylim([5000 9000])
title('Mean GMV by Site / After Combat')
imagesc(corr(gmv_data'));title('Intersubject Correlation Before Combat');figure;imagesc(corr(residuals));title('Structural Covariance After Combat');title('Intersubject Correlation After Combat');

imagesc(corr(gmv_data));title('Structural Covariance Before Combat');figure;imagesc(corr(residuals'));title('Structural Covariance After Combat')
%% Regional GMV difference
% idx = strcmp(diagnosis, 'AD') | strcmp(diagnosis, 'CN');
% for i=1:size(residuals,1)
%     [p,tbl,stats] = anova1((residuals(i,idx)),diagnosis(idx), 'off');
%     ps_gmv(i)=p;
%     fs_gmv(i)=tbl{2, 5};
%     stats_gmv{i}=stats;
% end


% Check the group-level differences

for i=1:size(residuals,1)
    [h,p,ci,stats] = ttest2((residuals(i,strcmp(diagnosis, 'SMC') )),(residuals(i,strcmp(diagnosis, 'CN') )) );
    d=meanEffectSize((residuals(i,strcmp(diagnosis, 'SMC') )),(residuals(i,strcmp(diagnosis, 'CN') )));
    ds_gmv(i)=d.Effect;
    ps_gmv(i)=p;
    ts_gmv(i)=stats.tstat;
end
corrected_ps=fdr_bh(ps_gmv);
results_to_display=ts_gmv(cortical_index).*corrected_ps(cortical_index);
figure;plot_cortical(parcel_to_surface(results_to_display,'aparc_fsa5'), 'cmap', 'RdBu_r', 'surface_name', 'fsa5', 'color_range', [-1*max(abs(results_to_display)) max(abs(results_to_display))])
results_to_display=ts_gmv(subcortical_index).*corrected_ps(subcortical_index);
figure;plot_subcortical(results_to_display,'ventricles','False','cmap', 'RdBu_r','color_range', [-1*max(abs(results_to_display)) max(abs(results_to_display))])

results_to_display=abs(ds_gmv(cortical_index));
figure;plot_cortical(parcel_to_surface(results_to_display,'aparc_fsa5'), 'cmap', 'Reds', 'surface_name', 'fsa5')
results_to_display=abs(ds_gmv(subcortical_index));
figure;plot_subcortical(results_to_display,'ventricles','False','cmap', 'Reds')


%% Correlate GM data with PET data

n_subjects = size(residuals, 2);
n_pet = size(petdata_mat, 1);
r_mat = NaN(n_subjects, n_pet);  % 1632 x 13

% Transpose residuals to be: subjects x 83
residuals_t = residuals';

% Compute subject-wise correlations for each PET map
for p = 1:n_pet
    pet_map = petdata_mat(p, :)';  % 83 x 1

    for s = 1:n_subjects
        subj_map = residuals_t(s, :)';  % 83 x 1
        r_mat(s, p) = corr(subj_map, pet_map, 'rows', 'pairwise');
    end
end

% short_names: 13x1 cell array of PET map labels
short_names = {'5HT1a','5HT1b','5HT2a','D1','D2','DAT','FDOPA','GABAa','Mu opioid','NAT','SERT','VAChT','mGluR5'}';

% Z-transform correlations
z_mat = atanh(r_mat);  % Fisher z-transform

% Define diagnosis order explicitly
group_order = {'CN', 'SMC', 'EMCI', 'LMCI', 'AD'};

n_groups = numel(group_order);
colors = lines(n_groups);

figure;
for p = 1:n_pet
    subplot(3, 5, p);  % For 13 PET maps
    hold on;

    % For ANOVA, collect all data and group labels for this PET map
    all_y = [];
    all_group_labels = [];

    for g = 1:n_groups
        group_name = group_order{g};
        group_idx = strcmp(diagnosis, group_name);
        y = z_mat(group_idx, p);  % Fisher z for group g and PET map p

        % Accumulate for ANOVA
        all_y = [all_y; y];
        all_group_labels = [all_group_labels; repmat({group_name}, sum(group_idx), 1)];

        % Plot violin or boxchart
        try
            violinplot(y, repmat({group_name}, sum(group_idx), 1), ...
                       'ShowMean', true, 'ViolinColor', colors(g,:));
        catch
            boxchart(repmat(g, sum(group_idx), 1), y, 'BoxFaceColor', colors(g,:));
        end
    end

    % Run ANOVA on all groups for this PET map
    [p_anova, tbl, stats] = anova1(all_y, all_group_labels, 'off');  % 'off' suppresses the figure
    F_score = cell2mat(tbl(2,5));  % Extract F-statistic from ANOVA table

    % Add F score to the title
    title(sprintf('%s\nF = %.2f, p = %.3g', short_names{p}, F_score, p_anova), 'Interpreter', 'none');

    xticks(1:n_groups);
    xticklabels(group_order);
    xtickangle(45);
    ylabel('Fisher z(r)');
    grid on;
end

sgtitle('Z-transformed Correlation of Subject GMV with PET Maps by Diagnosis');





% Plot with significance stars
n_subjects = size(residuals, 2);
n_pet = size(petdata_mat, 1);
r_mat = NaN(n_subjects, n_pet);  % 1632 x 13

% Transpose residuals to be: subjects x 83
residuals_t = residuals';

% Compute subject-wise correlations for each PET map
for p = 1:n_pet
    pet_map = petdata_mat(p, :)';  % 83 x 1

    for s = 1:n_subjects
        subj_map = residuals_t(s, :)';  % 83 x 1
        gm_z = (subj_map - mean(subj_map, 1)) ./ std(subj_map, 0, 1);      % [82 x 3227]
        pet_z = (pet_map - mean(pet_map, 2)) ./ std(pet_subset, 0, 2);  % [13 x 82]

        r_mat(s, p) = corr(subj_map, pet_map, 'rows', 'pairwise');
    end
end

% short_names: 13x1 cell array of PET map labels
short_names = {'5HT1a','5HT1b','5HT2a','D1','D2','DAT','FDOPA','GABAa','Mu opioid','NAT','SERT','VAChT','mGluR5'}';

% Z-transform correlations
z_mat = atanh(r_mat);  % Fisher z-transform

% Define diagnosis order explicitly
group_order = {'CN', 'SMC', 'EMCI', 'LMCI', 'AD'};
n_groups = numel(group_order);
colors = lines(n_groups);

figure;
for p = 1:n_pet
    subplot(3, 5, p);
    hold on;

    % Gather all data and groups for ANOVA
    all_y = [];
    all_group_labels = [];

    for g = 1:n_groups
        group_name = group_order{g};
        group_idx = strcmp(diagnosis, group_name);
        y = z_mat(group_idx, p);

        all_y = [all_y; y];
        all_group_labels = [all_group_labels; repmat({group_name}, sum(group_idx), 1)];

        % Plot violin or boxchart
        try
            violinplot(y, repmat({group_name}, sum(group_idx), 1), ...
                       'ShowMean', true, 'ViolinColor', colors(g,:));
        catch
            boxchart(repmat(g, sum(group_idx), 1), y, 'BoxFaceColor', colors(g,:));
        end
    end

    % Run ANOVA (suppress figure)
    [p_anova, tbl, stats] = anova1(all_y, all_group_labels, 'off');
    F_score = cell2mat(tbl(2,5));

    % Post-hoc multiple comparisons (Tukey's HSD)
    c = multcompare(stats, 'Display', 'off');

    % c columns: [group1 group2 lowerCI diff upperCI pvalue]
    sig_level = 0.05;
    sig_pairs = c(c(:,6) < sig_level, 1:2); % pairs with p < 0.05

    % Map group names to x positions for plotting stars
    group_x = 1:n_groups;

    % Draw significance stars between pairs
    y_max = max(all_y) + 0.1;  % starting height for stars
    star_gap = 0.05;  % vertical gap between stars if multiple pairs

    star_y = y_max;
    for k = 1:size(sig_pairs,1)
        g1 = sig_pairs(k,1);
        g2 = sig_pairs(k,2);
        x1 = group_x(g1);
        x2 = group_x(g2);

        % Draw line connecting the two groups
        plot([x1 x1 x2 x2], [star_y star_y+star_gap star_y+star_gap star_y], '-k', 'LineWidth', 1.5);

        % Put star centered above the line
        text(mean([x1 x2]), star_y+star_gap*1.1, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);

        star_y = star_y + star_gap * 1.5; % increase height for next star if multiple
    end

    % Title with F and p
    title(sprintf('%s\nF = %.2f, p = %.3g', short_names{p}, F_score, p_anova), 'Interpreter', 'none');

    xticks(1:n_groups);
    xticklabels(group_order);
    xtickangle(45);
    ylabel('Fisher z(r)');
    grid on;
end

sgtitle('Z-transformed Correlation of Subject GMV with PET Maps by Diagnosis');

%% Delta analysis
[n_regions, n_subjects] = size(residuals);
n_pet = size(petdata_mat, 1);

z_mat_loo = nan(n_subjects, n_pet, n_regions);  % Preallocate

% Loop over each region to leave it out to calculate the correlation values
% after removing each region
for r = 1:n_regions
    include_idx = setdiff(1:n_regions, r);  % leave one region out

    % Subset data excluding region r
    gm_subset = residuals(include_idx, :);         % [82 x 3227]
    pet_subset = petdata_mat(:, include_idx);      % [13 x 82]

    % Z-score (mean-centering and variance normalization)
    gm_z = (gm_subset - mean(gm_subset, 1)) ./ std(gm_subset, 0, 1);      % [82 x 3227]
    pet_z = (pet_subset - mean(pet_subset, 2)) ./ std(pet_subset, 0, 2);  % [13 x 82]

    % Compute correlation: dot product over regions (dim=1)
    % Transpose gm_z to [3227 x 82], so we can compute dot product with pet_z [13 x 82]
    corr_matrix = (gm_z' * pet_z') / (length(include_idx) - 1);  % [3227 x 13]

    % Fisher z-transform
    z_mat_loo(:, :, r) = atanh(corr_matrix);
end


for i=1:size(z_mat_loo,1)
    for j=1:size(z_mat_loo,2)
            difmat(i,j,:)=squeeze(z_mat_loo(i,j,:))-z_mat(i,j);
    end
end

% Do the group level analysis on delta
for j=1:size(z_mat_loo,2)
for i=1:size(z_mat_loo,3)
    [h,p,ci,stats] = ttest2(difmat(strcmp(diagnosis, 'AD'),j,i),difmat(strcmp(diagnosis, 'CN'),j,i) );
    d=meanEffectSize(difmat(strcmp(diagnosis, 'AD'),j,i),difmat(strcmp(diagnosis, 'CN'),j,i));
    ds_gmv(i)=d.Effect;
    ps_gmv(i)=p;
    ts_gmv(i)=stats.tstat;
end
corrected_ps=fdr_bh(ps_gmv);
results_to_display=ts_gmv(cortical_index).*corrected_ps(cortical_index);
figure;plot_cortical(parcel_to_surface(results_to_display,'aparc_fsa5'), 'cmap', 'RdBu_r', 'surface_name', 'fsa5', 'color_range', [-1*max(abs(results_to_display)) max(abs(results_to_display))])
title(sprintf('Delta difference ADxCN for %s', short_names{j}), 'Interpreter', 'none');

results_to_display=ts_gmv(subcortical_index).*corrected_ps(subcortical_index);
figure;plot_subcortical(results_to_display,'ventricles','False','cmap', 'RdBu_r','color_range', [-1*max(abs(results_to_display)) max(abs(results_to_display))])
title(sprintf('Delta difference ADxCN for %s', short_names{j}), 'Interpreter', 'none');

% results_to_display=abs(ds_gmv(cortical_index));
% figure;plot_cortical(parcel_to_surface(results_to_display,'aparc_fsa5'), 'cmap', 'Reds', 'surface_name', 'fsa5')
% results_to_display=abs(ds_gmv(subcortical_index));
% figure;plot_subcortical(results_to_display,'ventricles','False','cmap', 'Reds')
% title(sprintf('Delta difference between AD and CN for %s', short_names{j}), 'Interpreter', 'none');

end

% results_to_display=abs(ds_gmv(cortical_index));
% figure;plot_cortical(parcel_to_surface(results_to_display,'aparc_fsa5'), 'cmap', 'Reds', 'surface_name', 'fsa5',...
%     'color_range', [0.01 0.015])
% results_to_display=abs(ds_gmv(subcortical_index));
% figure;plot_subcortical(results_to_display,'ventricles','False','cmap', 'Reds')


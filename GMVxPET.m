%% The script below is used to calculate the correlation between GMV and PET maps across AD, bvFTD, and MCI groups, and the effect of individual ROIs on these correlations
%% The files in the github repository and the NIFTI toolbox are necessary for this script to run without errors
%% Cemal Koba, May 2024



clear

% Load subjects' metadata and mean PET data across 83 PET maps 
load('/home/koba/Desktop/kramer/subjects.mat')
petdata_full=readtable('/home/koba/Desktop/kramer/JuSpace-JuSpace_v1.5/JuSpace_v1.5/txtfiles/parcellated_petfiles.csv');

% Define the PET maps of interest 
pet_of_interest=[1 4 6 10 11 13 14 15 18 20 24 25 30];
petdata=petdata_full(pet_of_interest,:);
petdata_mat=table2array(petdata(:,2:end));

% Load DK40 parcellation and a sample PET map for getting the info of the
% nifti files 
dk40=load_nii('/home/koba/Desktop/kramer/JuSpace-JuSpace_v1.5/JuSpace_v1.5/atlas/atlas_desikan_resampled.nii');
template=load_nii('/home/koba/Desktop/kramer/JuSpace-JuSpace_v1.5/JuSpace_v1.5/PETatlas/5HT1b_P943_HC22.nii');


% Demographics 
[TABLE,CHI2,P,LABELS] = crosstab(nuisance(1:150,1),subjects(1:150,3)) % sex
sum(nuisance(:,1))/214;
sum(nuisance(strcmp(subjects(:,3),'AD'),1))/sum(strcmp(subjects(:,3),'AD'),1)
sum(nuisance(strcmp(subjects(:,3),'FTD'),1))/sum(strcmp(subjects(:,3),'FTD'),1)
sum(nuisance(strcmp(subjects(:,3),'MCI'),1))/sum(strcmp(subjects(:,3),'MCI'),1)



[P,ANOVATAB,STATS] =anova1(nuisance(:,2),subjects(:,3)) % age 
mean(nuisance(strcmp(subjects(:,3),'AD'),2))
std(nuisance(strcmp(subjects(:,3),'AD'),2))
mean(nuisance(strcmp(subjects(:,3),'FTD'),2))
std(nuisance(strcmp(subjects(:,3),'FTD'),2))
mean(nuisance(strcmp(subjects(:,3),'MCI'),2))
std(nuisance(strcmp(subjects(:,3),'MCI'),2))
mean(nuisance(:,2))
std(nuisance(:,2))

mdl2=anova1(nuisance(:,3),subjects(:,3)) % TIV
mean(nuisance(strcmp(subjects(:,3),'AD'),3))
std(nuisance(strcmp(subjects(:,3),'AD'),3))
mean(nuisance(strcmp(subjects(:,3),'FTD'),3))
std(nuisance(strcmp(subjects(:,3),'FTD'),3))
mean(nuisance(strcmp(subjects(:,3),'MCI'),3))
std(nuisance(strcmp(subjects(:,3),'MCI'),3))
mean(nuisance(:,3))
std(nuisance(:,3))



% Remove the nuisance variables from GM data
for i =1:size(fsdata,2)
    [B,BINT,R] = regress(fsdata(:,i),[ones(size(fsdata,1),1) nuisance]);
    fsdata_res(:,i) = R;
end
% residual check 
aa=mean(fsdata_res,2);
[H,P,CI,STATS] = ttest(aa)
histogram(aa)
mean(aa)
std(aa)


% Correlate GM data with PET data
for i=1:size(petdata,1)
    for j=1:length(fsdata)
        [rho,pval]=(corr(fsdata_res(j,:)',(petdata_mat(i,:))','Type','Spearman'));
        corrs(j,i)=rho;
        pvals(j,i)=pval;
    end
end


%mean correlations for total sample and across groups 
for i=1:size(petdata,1)
  totalmeans(i,1)=mean(atanh(corrs(:,i)));
  [H,P,CI,STATS] = ttest(atanh(corrs(:,i)));
  totalmeans(i,3)=P;
  totalmeans(i,2)=STATS.tstat;
    
  ftdmeans(i,1)=mean(atanh(corrs(strcmp(subjects(:,3),'FTD'),i)));
  [H,P,CI,STATS] = ttest(atanh(corrs(strcmp(subjects(:,3),'FTD'),i)));
  ftdmeans(i,3)=P;
  ftdmeans(i,2)=STATS.tstat;

  admeans(i,1)=mean(atanh(corrs(strcmp(subjects(:,3),'AD'),i)));
  [H,P,CI,STATS] = ttest(atanh(corrs(strcmp(subjects(:,3),'AD'),i)));
  admeans(i,3)=P;
  admeans(i,2)=STATS.tstat;

  mcimeans(i,1)=mean(atanh(corrs(strcmp(subjects(:,3),'MCI'),i)));
  [H,P,CI,STATS] = ttest(atanh(corrs(strcmp(subjects(:,3),'MCI'),i)));
  mcimeans(i,3)=P;
  mcimeans(i,2)=STATS.tstat;
end


% Check the group-level differences
for i=1:size(petdata,1)
    [p,tbl,stats] = anova1((atanh(corrs(:,i))),(subjects(:,3)), 'off');
    ps(i)=p;
    fs(i)=tbl{2, 5};
end
% Draw the graph
tiledlayout(5,3)
for i=1:size(petdata,1)
    % Left axes
    ax1 = nexttile;
    boxchart(ax1,(atanh(corrs(:,i))),'GroupByColor',categorical(cell2mat(subjects(:,2)))) % Use abs and corr > 0.1?
    if i == 1
        legend({'FTD','AD','MCI'}, 'FontSize',6);
    end
    ylabel('Correlation')
    file=(strrep(table2array(petdata(i,1)),'_',' '));
    file=strrep(file,'.nii','');
    xticklabels(cell(file))
end
sig_maps=fdr_bh(ps);
%% Find the most important ROIs

% Get the new correlation for each PET map, and subject after removing an ROI
for i=1:size(petdata,1)
    for j=1:length(fsdata)
        for k=1:size(fsdata,2)
            fs_temp=fsdata_res(j,:)';
            pet_temp=petdata_mat(i,:)';
            fs_temp(k)=[];
            pet_temp(k)=[];
            [rho,pval]=corr(fs_temp,pet_temp,'Type','Spearman');
            corrs_cook(j,i,k)=rho;
            pvals_cook(j,i,k)=pval;
        end
    end
end

% Calculate the delta -- difference in correlation
for k=1:size(petdata,1)
    corr_recalc=squeeze(corrs_cook(:,k,:));
    petvec=corrs(:,k);
    for i=1:size(fsdata,2)
        difmat(k,i,:)=atanh(corr_recalc(:,i))-atanh(petvec); %abs or not????
    end
end

%% Calculate the intra-group differences -- removed
% for i=1:size(petdata,1)
%     meanregions(i,:)=mean(abs(squeeze(difmat(i,:,:)))');
%     meanregions_ftd(i,:)=mean(abs(squeeze(difmat(i,:,strcmp(subjects(:,3),'FTD'))))');
%     meanregions_ad(i,:)=mean(abs(squeeze(difmat(i,:,strcmp(subjects(:,3),'AD'))))');
%     meanregions_mci(i,:)=mean(abs(squeeze(difmat(i,:,strcmp(subjects(:,3),'MCI'))))');  
% end
% 
% % Check the significance of these means 
% for i=1:size(petdata,1)
%     means_ftd=(squeeze(difmat(i,:,strcmp(subjects(:,3),'FTD'))));
%     means_ad=(squeeze(difmat(i,:,strcmp(subjects(:,3),'AD'))));
%     means_mci=(squeeze(difmat(i,:,strcmp(subjects(:,3),'MCI'))));
%     for j=1:size(petdata_mat,2)
%        [H,P,CI,STATS] = ttest(means_ftd(j,:));
%        intragroup_ftd_t(i,j)=STATS.tstat;
%        intragroup_ftd_p(i,j)=P;
% 
%        [H,P,CI,STATS] = ttest(means_ad(j,:));
%        intragroup_ad_t(i,j)=STATS.tstat;
%        intragroup_ad_p(i,j)=P;
% 
%        [H,P,CI,STATS] = ttest(means_mci(j,:));
%        intragroup_mci_t(i,j)=STATS.tstat;
%        intragroup_mci_p(i,j)=P;
%     end
% end
% for i=1:size(petdata,1);ps_corr_ftd(i,:)=fdr_bh(intragroup_ftd_p(i,:));end
% for i=1:size(petdata,1);ps_corr_ad(i,:)=fdr_bh(intragroup_ad_p(i,:));end
% for i=1:size(petdata,1);ps_corr_mci(i,:)=fdr_bh(intragroup_mci_p(i,:));end
% imagesc(ps_corr_ftd.*meanregions_ftd)
% imagesc(ps_corr_ftd(sig_maps,:).*meanregions_ftd(sig_maps,:))
% 
% 
% 
% % Save the mean intra-group differences as nifti - this part is manual 
% template=load_nii('/home/koba/Desktop/kramer/JuSpace-JuSpace_v1.5/JuSpace_v1.5/PETatlas/5HT1b_P943_HC22.nii');
% results=template;
% dk40_1d=reshape(dk40.img,[],1);
% results_to_map=mean(intragroup_ftd_t); %,sum(ps_corr_ftd) mean(ps_corr_ftd.*mean(intragroup_ftd_t))
% template_1d=zeros(size(dk40_1d));
% for i=1:83
%     template_1d(dk40_1d==i)=results_to_map(i);
% end
% results.img=reshape(template_1d,size(template.img));
% save_nii(results,'results/FTD_meandifference.nii.gz')

%% Calculate the inter-group differences for delta 
for i=1:size(petdata,1)
    means_ftd=(squeeze(difmat(i,:,strcmp(subjects(:,3),'FTD'))));
    means_ad=(squeeze(difmat(i,:,strcmp(subjects(:,3),'AD'))));
    means_mci=(squeeze(difmat(i,:,strcmp(subjects(:,3),'MCI'))));
    for j=1:size(petdata_mat,2)
       [H,P,CI,STATS] = ttest2(means_ftd(j,:),means_ad(j,:));
       intergroup_ftd_vs_ad_t(i,j)=STATS.tstat;
       intergroup_ftd_vs_ad_p(i,j)=P;

       [H,P,CI,STATS] = ttest2(means_ftd(j,:),means_mci(j,:));
       intergroup_ftd_vs_mci_t(i,j)=STATS.tstat;
       intergroup_ftd_vs_mci_p(i,j)=P;

       [H,P,CI,STATS] = ttest2(means_ad(j,:),means_mci(j,:));
       intergroup_ad_vs_mci_t(i,j)=STATS.tstat;
       intergroup_ad_vs_mci_p(i,j)=P;
    end
end
for i=1:size(petdata,1);ps_corr_ftdad(i,:)=fdr_bh(intergroup_ftd_vs_ad_p(i,:));end
for i=1:size(petdata,1);ps_corr_ftdmci(i,:)=fdr_bh(intergroup_ftd_vs_mci_p(i,:));end
for i=1:size(petdata,1);ps_corr_admci(i,:)=fdr_bh(intergroup_ad_vs_mci_p(i,:));end
sum(ps_corr_ftdad)
imagesc((ps_corr_ftdad(sig_maps,:).*mean(intergroup_ftd_vs_ad_t(sig_maps,:))))
figure
imagesc((ps_corr_ftdad.*mean(intergroup_ftd_vs_ad_t)))

mean(ps_corr_ftdmci)
mean(ps_corr_ftdmci.*mean(intergroup_ftd_vs_mci_t))

mean(ps_corr_ftdmci)
mean(ps_corr_ftdmci.*mean(intergroup_ftd_vs_mci_t))



%%  Brainstem volume difference between FTD and AD 

[H,P,CI,STATS] = ttest2(fsdata_res(strcmp(subjects(:,3),'FTD'),83),fsdata_res(strcmp(subjects(:,3),'AD'),83))

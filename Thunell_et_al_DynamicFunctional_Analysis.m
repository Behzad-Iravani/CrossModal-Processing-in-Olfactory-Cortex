% -*- UTF-8 -*-
% This script perform dynamic functional connnectivty  
% and network science analysis and is part of the analysis for
% "Whither unisensory olfactory cortex:
% processing of visual and auditory stimuli in olfactory cortex,
% independently of odor associations"  
%
% Copyright (C) Behzad Iravani
% behzadiravani@gmail.com
% 
% Department of Neurology and Neurological sciences, Stanford University, Palo Alto 
% 
% November, 2023 -- Philadelphia
% ------------------------------------------------------------------------

clear; clc; % clear memormy
% The preprocessing and denosing of the MRI data were carried out with CONN toolbox 22.a
% The denoised time-series from 272 ROIs were extracted for dynamic
% functional connectivity analysis 

% condition 0: intercept
% condition 1: pics
% condition 2: sounds

%-------------------------------------------------------------------------
clear all
clear  dFC

atlas_v = spm_vol('D:\SC09newData\atlas\P272.img');
[a.y, a.xyz] = spm_read_vols(atlas_v);
% retrive centroids
for j = 1:272
  ind = find(a.y == j);
  XYZ(:,j) = median(a.xyz(:, ind), 2); % centroids
  XYZa(:,j) = unique(a.y(a.y == j)); % centroids
end
%% estimate the dynamic functional connectivity
for isub = 1:47
    for icond = 1:2
        load(['data\' sprintf('ROI_Subject%03d_Condition%03d.mat', isub, icond)], 'data', 'names', 'conditionweights')
        FC = dyfc(cat(2,data{:,4:281}), names(:,4:281), conditionweights{1});
        FC = FC.connectivity();
       if icond == 1
           dFC.pic{isub} = FC;
       else
           dFC.snd{isub} = FC;
       end
    end % for icond
end % for isub
save results\dFC.mat dFC
%% load Gray Matter data
load gm.mat
gm = mean(cat(1,gm{:}));
%% TCs
plot(5*dFC.pic{7}.TC(:,1:10:end) + repmat(1:10:size(dFC.pic{7}.TC,2), size(dFC.pic{7}.TC,1), 1));
axis off
print Figures\resources\TCs.svg -vector -dsvg
% % ----
clf 
hold on
for steps = 100:110
    plot(floor(1+(steps-1)*10):floor(5 + (steps-1)*10), ...
        10*squeeze(dFC.pic{7}.wTC(:, steps, 1:10:end)) .* ...
        repmat(hann(5), 1, size(dFC.pic{7}.wTC(:, steps, 1:10:end),3)) ...
        + repmat(1:10:size(dFC.pic{7}.wTC,3), size(dFC.pic{7}.wTC,1), 1));
end

axis tight off
print Figures\resources\wTCs.svg -vector -dsvg
% % ----
for i=31:40
    subplot(1,10,i-30)
    imagesc(squeeze(dFC.pic{7}.dfc(:,:, i)))
    axis square
    colormap(viridis_white)
    axis off

end
print Figures\resources\dFC.svg -vector -dsvg
% %---


deletemepic = load(['data\' sprintf('ROI_Subject%03d_Condition%03d.mat', 7, 1)], 'data', 'names', 'conditionweights');
deletemesnd = load(['data\' sprintf('ROI_Subject%03d_Condition%03d.mat', 7, 2)], 'data', 'names', 'conditionweights');

imagesc([deletemepic.conditionweights{1}, deletemesnd.conditionweights{1}])
colormap(hot)
pbaspect([.3,1,1])
axis off
print Figures\resources\designMX.svg -vector -dsvg

%% second level statistics
pics = cellfun(@(x) x.beta, dFC.pic, 'UniformOutput', false);
snds = cellfun(@(x) x.beta, dFC.snd, 'UniformOutput', false);

pics = cat(3,pics{:});
snds = cat(3,snds{:});
tval.pics = (mean(pics,3)./(std(pics,[],3)/sqrt(size(pics,3))));
tval.snds = (mean(snds,3)./(std(snds,[],3)/sqrt(size(snds,3))));
% AVERAGE OVER HEMI

tval.pics = .5*(tval.pics(find([XYZ(1,:)>0,[0 0 0 1 1 1]]), find([XYZ(1,:)>0,[0 0 0 1 1 1]])) + ... % right 
tval.pics(find([XYZ(1,:)<0, [1 1 1 0 0 0]]),find([XYZ(1,:)<0, [1 1 1 0 0 0]]))); % left

tval.snds = .5*(tval.snds(find([XYZ(1,:)>0,[0 0 0 1 1 1]]), find([XYZ(1,:)>0,[0 0 0 1 1 1]])) + ... % right 
tval.snds(find([XYZ(1,:)<0, [1 1 1 0 0 0]]),find([XYZ(1,:)<0, [1 1 1 0 0 0]]))); % left
%% -------Visulization---------
subplot(121)
imagesc(tval.pics.*~eye(size(tval.pics)))
clim([0,30])
axis square off
title('\rm Pictures')
subplot(122)
imagesc(tval.snds.*~eye(size(tval.snds, 1,2)))
axis square off
title('\rm Sounds')
clim([0,30])
colormap("turbo")
axis off% sgtitle('\rm dFC t-values')
print Figures\resources\BetaMAP.svg -vector -dsvg

% remove those with low gray (40%) matter prob
tval.pics(gm<=.4, :) = 0;
tval.pics(:, gm<=.4) = 0;
tval.snds(gm<=.4,:) = 0;
tval.snds(:, gm<=.4) = 0;

subplot(121)
imagesc(tval.pics.*~eye(size(tval.pics, 1,2)))
clim([0,30])
axis square
title('\rm Pictures')
subplot(122)
imagesc(tval.snds.*~eye(size(tval.snds, 1,2)))
axis square
title('\rm Sounds')
clim([0,30])
sgtitle('\rm dFC group level t-stattistics')

print -dsvg -vector results\secondlevel.svg 
%% Network science analysis 
addpath('C:\MatlabToolboxes\BCT\2019_03_03_BCT')
tmpa = XYZa(:, XYZ(1,:)>0);
tmpx = XYZ(:, XYZ(1,:)>0);

A = threshold_proportional(tval.pics, .2);
[SPL,hops,Pmat] = distance_wei_floyd(A, 'inv');

retrieve_shortest_path(137,138, hops, Pmat) % LOC -> PPC 134 AMY
% retrieve_shortest_path(139,138, hops, Pmat) % AC -> PPC  67 eye movemnet


A = threshold_proportional(tval.snds, .2);
[SPL,hops,Pmat] = distance_wei_floyd(A, 'inv');

% retrieve_shortest_path(137,138, hops, Pmat) % LOC -> PPC 67 eye movement
retrieve_shortest_path(139,138, hops, Pmat) % AC -> PPC 134 AMY

%% result
out = double(a.y==tmpa(134));
outvol = atlas_v;
encode = 1;
for ROIs = ["rR_LOC", "rR_PPC", "rR_primaryAuditory"]
    encode =encode +1;
r_v = spm_vol(char(strcat("D:\SC09newData\atlas\ROIs\", ROIs, ".nii")));
[r.y, r.xyz] = spm_read_vols(r_v);
out(r.y>.9) = encode;
end

outvol.fname = 'results\output.nii';
spm_write_vol(outvol, out)



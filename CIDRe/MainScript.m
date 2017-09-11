%==========================================================================
%This MatLab-script shows usage and parameterization of the CIDRe-resampling
% and the min. Dist resampling of the synthtic example of paper GEO-2015-0220: 
% "CIDRe - a parameter constrained irregular resampling method for scattered point data"
%
% author:       Dipl.-Geoinf. Peter Menzel 
% institution:  Christian-Albrechts-Universität zu Kiel, Germany,
%               Department for Geosciences
%
% Example data: 
%  syntheticExample.csv:        contains the sythetic data set, used in
%                                   Section "Validation by synthetic data"
%   
% Additional MatLab-Functions:
%  CIDRe2D.m:                   contains CIDRe-resampling algorithm
%  Other_Resampling_Methods.m:  contains the additional resampling methods
%                                   used in the paper
%
% note: Please make sure that all needed files
%           - MainScript.m
%           - CIDRe2D.m:
%           - Other_Resampling_Methods.m
%           - syntheticExample.csv
%       are located in the same folder (cur. Matlab working directory). 
%       Execute script 
%           MainScript 
%       in the Matlab Command Window or open the file in the Matlab-Editor
%       and run (F5).
%
%==========================================================================
%
% By using this software, you agree to the terms of the SEG legal disclaimer, which may
% be found at: http://software.seg.org/disclaimer.txt
%
% This software may be found at: http://software.seg.org/2016/0002
%
%==========================================================================

%% 0) initialize, clear workspace, close figures etc.
clc;
clear all;
close all;
format long;
warning off;

%% 1) loading the data set
data_org = dlmread('syntheticExample.csv');
% prepare 3D coordinates - 3D coordinates are needed for CIDRe2D
data_org = [data_org(:,1:2) zeros(length(data_org),1) data_org(:,3)];

%% 2) Parametrizing and applying CIDRe
% Note 1: The results may differ slighly (not significant) from paper 
%           because of randomization in the kMeans-based resampling.
% Note 2: This way parameterizing CIDRe is only be used in this example ... 
%           currently a GUI-software (GUI-CIDRe) is used to define the parameters
%           and run the algorithm in practice.
% Note 3: The code is not fully optimzed and not parallized yet. 
%           For the used example a over all calculation time of <2000s 
%           will be needed on common hardware.
%
% parameters for result, shown in Fig. 5c
%for definition of these parameters please see header of CIDRe2D.m via
help CIDRe2D

RPW=  1.000000;
Wmethode= 1.000000;     
startRes= 80.000000;    
maxIter= 10.000000;
maxError= 10.000000;
minIRed= 0.950000;
keepBorder= 0.000000;
keepErrPts= 0.000000;
weightingOffset= 0.150000;
weightingScale= 1.000000;
ErrCalc= 0.000000;
ResamplingMethode= 2.000000;
% for parameterization of results for Fig. 4 adjust only parameters startRes and weightingOffset:
%startRes= 80;weightingOffset= 0.2;
%startRes= 80;weightingOffset= 0.15; %used example
%startRes= 80;weightingOffset= 0.1;
%startRes= 43;weightingOffset= 0.1;
%startRes= 43;weightingOffset= 0.05;
%startRes= 80;weightingOffset= 0.0;
%startRes= 31;weightingOffset= 0.0;
%startRes= 16;weightingOffset= 0.0;
%startRes= 10;weightingOffset= 0.0;
%startRes= 7;weightingOffset= 0.0;
%startRes= 7;weightingOffset= -0.05;

CIDRe_result = CIDRe2D(data_org,...
                     RPW,...                        
                     Wmethode,...
                     startRes,...
                     maxIter,...
                     maxError,...
                     minIRed,...
                     keepBorder,...
                     keepErrPts,...
                     true,...                       %please set to false if no progress output is wanted
                     weightingOffset,...
                     weightingScale,...
                     ErrCalc,...
                     ResamplingMethode,...
                     1);                            %always use this value

                 
%% 3) Parametrizing and applying minimum distance resampling 
% parameters for result, shown in Fig. 5a
stepSize = 2;
% for parameterization of results for Fig. 4
%stepSize = 0.5;
%stepSize = 1;
%stepSize = 1.5;
%stepSize = 2; %used example
%stepSize = 2.5;
%stepSize = 5;
%stepSize = 10;
minDist_result = Other_Resampling_Methods(data_org(:,1:2),data_org(:,4),2,stepSize,'minDist',stepSize); 

%% 4) Restoring the original data set and calculating the RMSE
% using Matlab Natural-Neigthbor-Interpolation

% restoring the CIDRe results, Fig. 5c
interpa = TriScatteredInterp(CIDRe_result(:,1),CIDRe_result(:,2),CIDRe_result(:,4),'natural');
tmp = interpa(data_org(:,1),data_org(:,2));
% removing extrapolation-caused NaN by Nearest-Neightbor-Results
interpb = TriScatteredInterp(CIDRe_result(:,1),CIDRe_result(:,2),CIDRe_result(:,4),'nearest');
tmp(~isfinite(tmp)) = interpb(data_org(~isfinite(tmp),1),data_org(~isfinite(tmp),2));
CIDRe_restored = tmp;   
% RMSE
CIDRe_diff = data_org(:,4)-CIDRe_restored;
CIDRe_RMSE = sqrt(mean(CIDRe_diff.*CIDRe_diff));
display(['RMSE for CIDRe: ',num2str(CIDRe_RMSE)]);


% restoring the CIDRe results, Fig. 5b
interpa = TriScatteredInterp(minDist_result(:,1),minDist_result(:,2),minDist_result(:,3),'natural');
tmp = interpa(data_org(:,1),data_org(:,2));
% removing extrapolation-caused NaN by Nearest-Neightbor-Results
interpb = TriScatteredInterp(minDist_result(:,1),minDist_result(:,2),minDist_result(:,3),'nearest');
tmp(~isfinite(tmp)) = interpb(data_org(~isfinite(tmp),1),data_org(~isfinite(tmp),2));
minDist_restored = tmp;   

% RMSE
minDist_diff = data_org(:,4)-minDist_restored;
minDist_RMSE = sqrt(mean(minDist_diff.*minDist_diff));
display(['RMSE for min. Dist: ',num2str(minDist_RMSE)]);

%% 5) plot the data to reproduce figures 5 a,b,c,d
%create colormap used in the paper
tmp =([0:0.01:0.9].^0.5)';
errCM = [tmp tmp tmp;repmat([0.95 0.95 0.95],23,1);tmp(end:-1:1) tmp(end:-1:1) tmp(end:-1:1)];

%get the figure limits
loc = data_org(:,2) >=200 & data_org(:,2) <=400;

%5a) min. Dist. result
locR1 = minDist_result(:,2) >=200 & minDist_result(:,2) <=400;
figure('name','a) min. Dist. result','OuterPosition',[0 0 1600 500],'Position',[0 0 1600 500],'MenuBar','none','color',[1 1 1]);
set(gca,'position',[0 0 0.8 1])
scatter(minDist_result(locR1,1),minDist_result(locR1,2),2,minDist_result(locR1,3),'filled')
daspect([1 1 1])
axis off
caxis([0 300])
colormap(gray)
cb=colorbar('fontsize',20,'YTick',[50 150 250]);
ylabel(cb,'synthetic parameter','fontsize',20)
xlim([0 500])
ylim([200 400])

%5b) min. Dist. differences
figure('name','b) min. Dist. difference','OuterPosition',[0 0 1600 500],'Position',[0 0 1600 500],'MenuBar','none','color',[1 1 1]);
set(gca,'position',[0 0 0.8 1])
scatter(data_org(loc,1),data_org(loc,2),2,minDist_diff(loc,1),'filled')
daspect([1 1 1])
axis off
caxis([-20 20])
colormap(errCM)
cb=colorbar('fontsize',20,'YTick',[-20 -10 0 10 20]);
ylabel(cb,'difference','fontsize',20)
xlim([0 500])
ylim([200 400])

%5c) CIDRe result
locR2 = CIDRe_result(:,2) >=200 & CIDRe_result(:,2) <=400;
figure('name','c) CIDRe result','OuterPosition',[0 0 1600 500],'Position',[0 0 1600 500],'MenuBar','none','color',[1 1 1]);
set(gca,'position',[0 0 0.8 1])
scatter(CIDRe_result(locR2,1),CIDRe_result(locR2,2),2,CIDRe_result(locR2,4),'filled')
daspect([1 1 1])
axis off
caxis([0 300])
colormap(gray)
cb=colorbar('fontsize',20,'YTick',[50 150 250]);
ylabel(cb,'synthetic parameter','fontsize',20)
xlim([0 500])
ylim([200 400])

%5d) CIDRe differences
figure('name','d) CIDRe difference','OuterPosition',[0 0 1600 500],'Position',[0 0 1600 500],'MenuBar','none','color',[1 1 1]);
set(gca,'position',[0 0 0.8 1])
scatter(data_org(loc,1),data_org(loc,2),2,CIDRe_diff(loc,1),'filled')
daspect([1 1 1])
axis off
caxis([-20 20])
colormap(errCM)
cb=colorbar('fontsize',20,'YTick',[-20 -10 0 10 20]);
ylabel(cb,'difference','fontsize',20)
xlim([0 500])
ylim([200 400])

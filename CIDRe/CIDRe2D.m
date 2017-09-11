function [reducedData,AbwPPt,iData] = CIDRe2D(origData,RPW,Wmethode,startRes,maxIterations,maxError,maxReductFactor,doKeepBorder,doKeepErrPts,outputProgress,weightingOffset,weightingScale,doErrCalc,resMeth,olm)
%=========================================================================
%(C)ontrained (I)ndicator (D)ata (Re)sampling for horizontal (2d) PointData in 3d-Space
%This function contains the algorithm for CIDRe, which was used to resampled the examples in paper GEO-2015-0220: 
% "CIDRe - a parameter constrained irregular resampling method for scattered point data"
%%Revised version of the code.
%
%   author:       Dipl.-Geoinf. Peter Menzel
%   institution:  Christian-Albrechts-Universitaet zu Kiel, Germany, Department for Geosciences
%   date:         2015-04 (revision 2015-10)
%
% By using this software, you agree to the terms of the SEG legal disclaimer, which may
% be found at: http://software.seg.org/disclaimer.txt
%
% This software may be found at: http://software.seg.org/2016/0002
%
%--------------------------------
%
% [reducedData,AbwPPt,iData] = CIDRe2D(origData,...
%                                RPW,Wmethode,startRes,maxIterations,maxError,maxReductFactor,...
%                                doKeepBorder,doKeepErrPts,outputProgress,weightingOffset,...
%                                weightingScale,doErrCalc,resMeth,olm)
%   input:
%       origData        = input DataPoints [x y z (param1 ... paramN)]
%       RPW             = global relativ Parameterweights [wP1 ... wPN],
%                           only used for data sets with >1 parameter 
%                           values per point
%       Wmethode        = 1 or default=strategy I, 2=strategy II, the
%                           proposed strategies Ia and IIb are not included
%                           in this version
%       startRes        = start parameter for sectorization 
%       maxIterations   = max. number of iterations (or < 0: interation count is not abort-criterium)
%       maxError        = Max. error limit (%)
%       maxReductFactor = max. reduction factor, if inter-interation-
%                           reduction rate is higher->abort
%       doKeepBorder    = keep bordering points
%       doKeepErrPts    = keep point above error limit
%       outputProgress  = print progress
%       weightingOffset = global offset for limiting weights to control
%                           resampling results, as higher value, as higher
%                           resampling rate (corresponds to value W^{offset} for calculation of the indicator weights)
%       weightingScale  = scale the calculated weights to control the "shape"
%                           of the histogram ... needed for "bad"
%                           weight-distribution
%       doErrCalc       = perform error calculation (not needed for the example)
%       resMeth         = resampling method, always used as 2 (irregular
%                           kMeams-based), but other methods are provided
%                           as well
%       olm             = "overlaymethod", only needed for multi-parameter weighting, always used as 1!
%   output:
%       reducedData     = reduced dataset
%       AbwPPt          = difference for each dismissed point (0 for used
%                           points) -> empty if doErrCalc==false
%       iData           = interpolated(restored) data -> empty if doErrCalc==false
%
%Notes:
%1) This version of CIDRe does NOT include some features, proposed in the paper:
%   -weighting strategy Ia and IIa
%   -outlier detection
%   -3D point cloud resampling
% These features are not used for the examples, shown in the Sections "Validation by synthetic data"
%   and "Real data application: ..." and will be integrated in future 
%   versions of CIDRe.
%
%2) Multi-parameter resampling is at least possible, but not fully tested
%       and may produce "not intuitive" results for real data.
%       This feature is currently under developement and not needed/used
%       for the examples, as well.
%=========================================================================

%% 0) initial step
display('CIDRe2D for paper GEO-2015-0220 ')
recExt = 1;
ErrEstM = 3;

display(' ')
display('CIDRe2D(...) ... begin data resampling ...')
refDim = 4;%column 4 and higher are parameter values
dataPoints = origData;
s = size(dataPoints);
sOriginalData = size(origData);
origParaRange  = 1;

abortByIter = false;
if maxIterations > 0
    abortByIter = true;
end

%RPW will be always 1 for the examples, it is only needed for mutli-parameters
if size(origData,2)-3 < size(RPW,2)
    RPW = RPW(1:size(origData,2)-3);
else
    RPW = [RPW ones(1,(size(origData,2)-3)-size(RPW,2))];
end
%normalize paramter values
pMax = max(origData(:,4:end));
pMin = min(origData(:,4:end));
pRange = pMax-pMin;
nWParam = size(origData,2)-3; %always 1 for the examples used in the paper
for np = 1:nWParam 
    origData(:,3+np)=((origData(:,3+np)-pMin(np))./pRange(np));
end
overlayMethode =olm;%1 or 0; 1 = superposition, 0 = maximum
dataPoints = origData;%use normalized data

iteration = 1;
endIt = maxIterations;
useMapParam = 1;

AbwPPt = zeros(s(1),nWParam);
RelErrPPt = zeros(s(1),nWParam);

sumError = 0;
ErrorToBig = false;
noSigDiff = false;
red = 0;
maxKeepW = -1;

minC = startRes;
if minC < 5
    minC = 5;
end
display(['CIDRe2D(...) .... minimum sector containment: ',num2str(minC)]);
%% iteration loop
while not(ErrorToBig) && not(noSigDiff) && (not(abortByIter) || iteration <= endIt)
    display(['CIDRe2D(...) ... ',num2str(iteration),'. iteration']);
    %% 1) preprocessing step

    display('CIDRe2D(...) .... search-preprocessing ....')
    resP = zeros(length(dataPoints),1);
    tic
    %compute the discretisation for sectoring
    sP = size(dataPoints);
    ratio = (max(dataPoints(:,1))-min(dataPoints(:,1)))/(max(dataPoints(:,2))-min(dataPoints(:,2)));
    %sector resolution:
    sx = ceil(sqrt(ratio*sP(1)/minC));
    sy = ceil(sx/ratio);
    %create sectoring and search-neightborhood
    [searchIDs,nMap,ss,ns,origin] = createNeighborhoodMapAdaptive_local(dataPoints(:,1:2),[sx sy],1);
    %searchIDs contains the sector-id for each data point ... is not needed
    %   for the following steps.
    %nMap containes for each sector all id's of neightboring points
    %ss,ns,origin contain sectoring informations
    %ss corresponds to the sector size, introduced in Section 'CIDRe – Constrained Indicator Data Resampling', Subsection 'Preprocessing'
    sx = ns(1);
    sy = ns(2);
    display(['CIDRe2D(...) .... search-origin: ',num2str(origin)])
    display(['CIDRe2D(...) .... search-stepsize: ',num2str(ss)])    
    display(['CIDRe2D(...) .... search-sectors: ',num2str(ns)])
    r = mean(ss);
    %r corresponds to R_max, introduced in Section 'CIDRe – Constrained Indicator Data Resampling', Subsection 'Preprocessing'
    toc    
    %% 2) weighting step
    display('CIDRe2D(...) .... calculate parameter weights ....') 
    tic
    NLP = cell2mat(nMap(:,1));
    idim = size(NLP);    
        %loop for all sectors
        for i = 1:idim(1)
            pointIds = cell2mat(nMap(i,3));
            sektorPoints = dataPoints(pointIds,:);
            npointIds = cell2mat(nMap(i,2));
            neigthborPoints = dataPoints(npointIds,:);
            
            ssp = size(sektorPoints);
            snp = size(neigthborPoints);
            %loop for all sector points
            for sp = 1:ssp(1)        
                p = sektorPoints(sp,:);
                pID = pointIds(sp);
                nCount = 0;
                %loop for all neighboring points of a sector
                for np = 1:snp(1)
                    q = neigthborPoints(np,:);
                    dist = norm(p(1:3)-q(1:3));
                    if dist > 0 && dist <= r
                        nCount = nCount + 1;                     
                        
                        %loop for all parameter values (1 for given example)
                        pW = zeros(1,nWParam);
                        for usedP = 1:nWParam
                            pGrad4 = (abs((p(usedP+3)-q(usedP+3)))/dist); %dp_N for each neighbor, Equation 1
                            pW(usedP) = pGrad4*RPW(usedP);
                        end
                        tmp = 0;
                        for usedP = 1:nWParam                                
                            if overlayMethode == 0
                                if tmp < pW(usedP)
                                    tmp = pW(usedP);
                                end
                            else
                               tmp = tmp + pW(usedP);
                            end
                        end 

                        if ~(Wmethode == 2)%alsways use arithm. mean except for Wmethode == 2
                            %"arithm. mean"-weighting - Strategy I, Equation 2
                            resP(pID,useMapParam) = resP(pID,useMapParam)+tmp;
                        else
                            %RMS-weighting - Strategy II, Equation 3
                            resP(pID,useMapParam) = resP(pID,useMapParam)+(tmp*tmp); % Summize the squared value
                        end
                    end            
                end
                if nCount > 0
                    %calculate the mean, Equation 2 and Equation 3
                    resP(pID,useMapParam) = resP(pID,useMapParam)/(nCount*(nWParam^overlayMethode));
                    %root the squared values
                    if Wmethode == 2
                        %only needed for RMS, Equation 3
                        resP(pID,useMapParam) = sqrt(resP(pID,useMapParam));
                    end
                end
            end
            
            if outputProgress
                printProgress(toc,i,idim(1));
            end
        end
        toc
        %scale weights by iteration
        resP = resP.^(1/iteration); 
        %normalized weights
        resP(:,useMapParam) = (resP(:,useMapParam)-min(resP(:,useMapParam)))/(max(resP(:,useMapParam))-min(resP(:,useMapParam))); 
        %scale the weights weights
        if weightingScale >=0.2 && weightingScale <=1.5
            resP(:,useMapParam) = resP(:,useMapParam)/weightingScale;
            resP((resP(:,useMapParam)>1),useMapParam) =1;
        end
                     
        %% calculate the initial limiting weight based on the 95% quantil
        histq = 100;
        if iteration == 1
            h = hist(resP(:,useMapParam),histq);
            N = length(resP(:,useMapParam));         
            n = 1;
            while sum(h(1:n)) < 0.95*N
                n= n+1;
            end
            maxKeepW = n/length(h);
            if maxKeepW > 0.9
                maxKeepW = 0.9;
            end
        end       
        %maxKeepWeight corresponds to W^{95%}, introduced in Section 'CIDRe – Constrained Indicator Data Resampling', Subsection 'Irregular resampling'
    %% 3) irregular resampling based on the caluclated point weights
    display('CIDRe2D(...) ... point resampling ...')   
    tic
    newPoints = [];
    simpleStatus = zeros(sP(1),1);
    seks = cell2mat(nMap(:,1));

    if doKeepBorder
        isBS = getBorderSektors_local(seks);
    else
        isBS = ones(size(seks,1),1) == 0;
    end
        
    maxOff = weightingOffset;
    %maxOff corresponds to W^{offset}, introduced in Section 'CIDRe – Constrained Indicator Data Resampling', Subsection 'Irregular resampling'
    wOffset = maxOff-0.05*(iteration-1);

    display(['CIDRe2D(...) .... Resampling with limiting weight ',num2str(maxKeepW+wOffset)]);
    
    newPoints = zeros(size(dataPoints));
    nnPts = 0;
    %loop for all sectors - resample each sector separatly
    for i = 1:idim(1)
        doRed = true;
        if doKeepBorder && isBS(i) == 1
            doRed = false;
        end
        pointIds = cell2mat(nMap(i,3));
        sektorPoints = dataPoints(pointIds,:);
        
        if doRed
            bVal = maxKeepW+wOffset;
            if bVal > 0.9
                bVal = 0.9;
            elseif bVal < 0.1
                bVal = 0.1;
            end
            %bval corresponds to W_{I}^{ref}, introduced in Section 'CIDRe – Constrained Indicator Data Resampling', Subsection 'Irregular resampling'
            %select the used points from the given data points
            addPoints = resampleByWeighting(sektorPoints,resP(pointIds,useMapParam),bVal,4,2,resMeth);
            %remove doubled points - should not be required, but may
            %prevent problems
            addPoints = unique(addPoints,'rows');
            newPoints(nnPts+1:nnPts+size(addPoints,1),:) = addPoints;
            nnPts = nnPts+size(addPoints,1);
        else
            newPoints(nnPts+1:nnPts+size(sektorPoints,1),:) = sektorPoints;
            nnPts = nnPts+size(sektorPoints,1);
        end     
         
        if outputProgress
            printProgress(toc,i,idim(1));
        end    
    end
    newPoints = newPoints(1:nnPts,:);
    toc
   
    %%  iteration check, is a next iteration possible/needed?   
    Reduktionsfaktor = length(newPoints)/length(dataPoints);
    display(['CIDRe2D(...) .... reductionsfactor: ',num2str(Reduktionsfaktor),', Iteration: ',num2str(iteration)])
    if Reduktionsfaktor > maxReductFactor
            noSigDiff = true;
            display('CIDRe2D(...) .... finished, no significant changes!')
    else
        if not(abortByIter) || iteration + 1 <= endIt
            dataPoints = newPoints;
        else
            display('CIDRe2D(...) .... finished, iteration limit reached!')
        end
    end
    iteration = iteration +1;
end

if doErrCalc
    %% error estimation - not reliable and unefficiant in this version, use external calculations (see MainScript.m Part 4)
    %This part is not needed, because of doErrCalc=false for all examples in
    %the paper. The errors are calculated separately.    
    % interpolate to original locations and calculate differences
    %and find points with to large errors. These points shall be kept in
    %the resampled data if doKeepErrPts==true.
    r = recExt*r;
    display('CIDRe2D(...) ... analyze Errors ...')   
    
    %use previousely defined search sektors for finding neighboring points
    %to interpolate on the original point locations
     [origSearchIDs,origNMap,origSS,origSOrigin] = createNeighborhoodMap_local(origData(:,1:2),[sx sy],0);
     [newSearchIDs,newNMap] = createNeighborhoodMapByGivenSektorization_local(newPoints(:,1:2),origSOrigin,origSS,1);      

    tic
    origIDXs = cell2mat(origNMap(:,1));
    origIDim = size(origIDXs);
    newIDXs = cell2mat(newNMap(:,1));
    newIDim = size(newIDXs);
    interpProp = zeros(sOriginalData(1),nWParam); %restored parameter values
    noInterpPts = false(size(origData,1),1);
    RMS = zeros(sOriginalData(1),1); %root-mean-squered-Error
    %interpolated resampled data to old locations
    for i = 1:origIDim(1)
        opIds = cell2mat(origNMap(i,2));
        origPoints = origData(opIds,:);

        nNList = find(ismember(newIDXs,origIDXs(i,:),'rows'));
        if not(isempty(nNList))
            reducedPoints = newPoints(cell2mat(newNMap(nNList,2)),:);
            
            sRP = size(reducedPoints); 
            if sRP(1) > 2
                reducedTris = delaunay(reducedPoints(:,1),reducedPoints(:,2));

                for oP = 1:length(opIds)
                    matchFound = false;
                    nPidx = 1;
                    usedVals = zeros(sRP(1),1); %% weightings for IDW
                    nVals = 0;
                    sumW = 0;
                    mindist = r;
                    matchVal = 0;
                    if  ErrEstM == 1 || ErrEstM == 2   
                        while not(matchFound) && nPidx <= sRP(1)
                            dist = norm(origPoints(oP,1:3)-reducedPoints(nPidx,1:3));
                            distEps = 0;%0.1*r;
                            if dist > distEps && dist < r
                                weight = ((r-dist)/(r*dist))^2;
                                interpProp(opIds(oP)) = interpProp(opIds(oP)) + (weight*reducedPoints(nPidx,refDim));
                                sumW = sumW + weight;
                                usedVals(nPidx) = reducedPoints(nPidx,refDim);
                                if dist < mindist
                                    mindist = dist;
                                    matchVal = reducedPoints(nPidx,refDim);
                                end
                                nVals = nVals+1;
                                nPidx = nPidx+1;
                            elseif dist <= distEps
                               interpProp(opIds(oP)) = origPoints(oP,refDim);
                               matchFound = true;
                            else
                                nPidx = nPidx+1;
                            end
                        end
                        if not(matchFound)
                            if not(sumW == 0) && not(ErrEstM == 2)
                                interpProp(opIds(oP)) = interpProp(opIds(oP))/sumW;%%r-sphere IDW
                            else
                                interpProp(opIds(oP)) = matchVal;%%nearest-neightbor
                            end
                        end
                    elseif ErrEstM == 3
                        interpProp(opIds(oP),:) = interpolateNaN2D(origPoints(oP,1:3), reducedPoints, reducedTris);
                    end
                end
            else
                display(['CIDRe2D(...) .... Warning, to less points for good estimation! (',num2str(sRP(1)),')'])
                for oP = 1:length(opIds)
                    interpProp(opIds(oP),:) = mean(reducedPoints(:,4:end));
                end 
            end
        else
            display(['CIDRe2D(...) .... Warning, no points for estimation!'])
            for oP = 1:length(opIds)
                interpProp(opIds(oP),:) = min(origPoints(:,4:end)) + ((max(origPoints(:,4:end))-min(origPoints(:,4:end)))*0.5);                
            end   

        end
        if outputProgress
            printProgress(toc,i,origIDim(1));
        end
    end

    for usedP = 1:nWParam
        %globalize paramter values        
        origData(:,usedP+3) = origData(:,usedP+3)*pRange(usedP) + pMin(usedP);
        interpProp(:,usedP) = interpProp(:,usedP)*pRange(usedP) + pMin(usedP);
        newPoints(:,usedP+3) = newPoints(:,usedP+3)*pRange(usedP) + pMin(usedP);
    end 

    iData = interpProp;
    diff = origData(:,4:end)-interpProp(:,1:end); 
    
    toc        
    %keep points with large differences
    AbsAbw = zeros(size(diff));    
    if origParaRange > 0
        for usedP = 1:nWParam
            AbsAbw(:,usedP) = (abs(diff(:,usedP))/pRange(usedP))*100; 
        end

        if doKeepErrPts
            idxs = logical(zeros(size(AbsAbw,1),1));
            for usedP = 1:nWParam
                if RPW(usedP) > 0
                    loc = AbsAbw(:,usedP)>maxError;
                    idxs(loc) = true;
                end
            end

            errPts = origData(idxs,:);            
            errPtsS = size(errPts);            
            display(['CIDRe2D(...) .... ',num2str(errPtsS(1)),' points arent resampled because of large errors! ...'])            
            newPoints = [newPoints;errPts];
            
            AbsAbw(idxs,:) = 0;
            diff(idxs,:) = 0;
            iData(idxs,:) = origData(idxs,4:end);
        end  
    else
        AbsAbw = diff;
    end 

    sumError = zeros(1,size(AbsAbw,2));
    for usedP = 1:nWParam
        %arithmetic mean
        %sumError(1,usedP) = mean(AbsAbw(not(AbsAbw(:,usedP) == 0),usedP));
        %RMS
        sumError(1,usedP) = sqrt(mean(AbsAbw(not(AbsAbw(:,usedP) == 0),usedP).^2));
    end
        drchschnAbw = sumError;
        maxAbw = max(AbsAbw);
            display(['CIDRe2D(...) .... RMS rel. difference: ',num2str(drchschnAbw)])
            display(['CIDRe2D(...) .... max. rel. difference: ',num2str(maxAbw)])
                dataPoints = newPoints;        
                AbwPPt = diff;
            cummulatedReduction = (1-length(dataPoints)/sOriginalData(1))*100;
            display(['CIDRe2D(...) .... Recuction ',num2str(cummulatedReduction),'% -> to ',num2str(length(dataPoints)),' points'])
    display('CIDRe2D(...) ... ready ...')
else
    %% return the resampled data without error calculation
    for usedP = 1:nWParam
        %globalize parameter values
        origData(:,usedP+3) = origData(:,usedP+3)*pRange(usedP) + pMin(usedP);
        newPoints(:,usedP+3) = newPoints(:,usedP+3)*pRange(usedP) + pMin(usedP);
    end       
    dataPoints = newPoints; %returned result
    AbwPPt = [];%not calculated
    iData = [];%not calculated
    
    cummulatedReduction = (1-length(dataPoints)/sOriginalData(1))*100;
    display(['CIDRe2D(...) .... Reduction ',num2str(cummulatedReduction),'% -> to ',num2str(length(dataPoints)),' points'])
end
display(' ')
reducedData = dataPoints;
end
%% local functions for CIDRe2D
%==========================================================================
%bi-linear interpolation
function [values,noInterp] = interpolateNaN2D(iPoint, dataPoints, dataTris)
%bi-linear triangultion based interpolation using barycentric coordinates
sP = size(dataPoints);
sT = size(dataTris);
[ismem,idx] = ismember(iPoint,dataPoints(:,1:3),'rows');
noInterp = false;
if ismem
    %this point is contained in the resampled data set, no interpolation
    %needed
    values = dataPoints(idx,4:end);
    noInterp = true;
else
    found = false;
    t = 1;
    while not(found) && t <= sT(1)
        P = [iPoint(1:2) 0];
        A = [dataPoints(dataTris(t,1),1:2) 0];
        B = [dataPoints(dataTris(t,2),1:2) 0];
        C = [dataPoints(dataTris(t,3),1:2) 0];
                
        if(PinABC(P,A,B,C))            
            PA = P-A;
            BA = B-A;
            CA = C-A;
            
            %if norm(cross(BA,CA)) > 0.0001
            h1 = (PA(1)*BA(2)-PA(2)*BA(1))/((CA(1)*BA(2)-CA(2)*BA(1)));
            h2 = (PA(1)*CA(2)-PA(2)*CA(1))/((CA(2)*BA(1)-CA(1)*BA(2)));
            if isfinite(h1) && isfinite(h2)
                values = h1*dataPoints(dataTris(t,3),4:end) + h2*dataPoints(dataTris(t,2),4:end) + (1-h1-h2)*dataPoints(dataTris(t,1),4:end);
            else
            %arithmetic mean
                values=(1/3)*(dataPoints(dataTris(t,1),4:end)+dataPoints(dataTris(t,2),4:end)+dataPoints(dataTris(t,3),4:end));
            end   
            found = true;
        else
            t = t +1;
        end
    end
    %no containing triangle found, use nearest neightbor
    if not(found)
        mindist = norm(iPoint-dataPoints(1,1:3));
        id = 1;
        for i = 1: sP(1)
            dist = norm(iPoint-dataPoints(i,1:3));
            if dist < mindist
                mindist = dist;
                id = i;
            end
        end
        values = dataPoints(id,4:end);
    end
end
end
%point in triangle test
function [in] = PinABC( P,A,B,C)
% 3d-Point P in Triangle of ABC?
in = true;

% %3)own methode also implemented in c++, using max/min for scalars (6x faster then 2)
if P(1) > max(max(A(1),B(1)),C(1)) || P(1) < min(min(A(1),B(1)),C(1))
    in = false;
    %display('PinABC: point out on x-Axis');
elseif P(2) > max(max(A(2),B(2)),C(2)) || P(2) < min(min(A(2),B(2)),C(2))
    in = false;
    %display('PinABC: point out on y-Axis');
elseif P(3) > max(max(A(3),B(3)),C(3)) || P(3) < min(min(A(3),B(3)),C(3))
    in = false;        
    %display('PinABC: point out on z-Axis')
end

%second try 3 same-side-tests
if in == true
    if (not(IsSameSide(P,A,B,C)) || ...
            not(IsSameSide(P,B,A,C)) || ...
            not(IsSameSide(P,C,A,B)))
        in = false;
        %display('PinABC: point out cause of point-side-testing')
    else
        %display('PinABC: point inside or on one edge')
    end
end

end
%is same side test
function [ is ] = IsSameSide( P1, P2, A, B )
 BA = (B-A)/norm(B-A);
 P1A = (P1-A)/norm(P1-A);
 P2A = (P2-A)/norm(P2-A);
 
 v1 = cross(BA,P1A);
 v2 = cross(BA,P2A);
 
 if (v1(1) == 0 && v1(2) == 0 && v1(3) == 0) || (v2(1) == 0 && v2(2) == 0 && v2(3) == 0) 
     %one point on edge, true
     is = true;
 else
    v1 = v1/norm(v1);
    v2 = v2/norm(v2);
    DOT = dot(v1,v2);
    if DOT >= 0
        is = true;
    else
        is = false;
    end
 end

end
%resampling by point weights, 10% steps
function [resampledData] = resampleByWeighting(dataPoints, pWeights, maxW, varargin)
minRes = 4;
if nargin > 3
    minRes = varargin{1};
end
dim = 2;
if nargin > 4
    dim = varargin{2};
end
meth = 1;
if nargin > 5
    meth = varargin{3};
end

stepVals = [0.3*maxW 0.5*maxW maxW 1.2*maxW];
%stepVals corresponds to indicator weights W_{I}, described in Section 'CIDRe – Constrained Indicator Data Resampling', Subsection 'Irregular resampling'
redRate = [0.01 0.1 0.33 0.8];
%redRate corresponds to indicator reduction factors RF_{I}, described in Section 'CIDRe – Constrained Indicator Data Resampling', Subsection 'Irregular resampling'

undoneData = dataPoints;
undoneWeights = pWeights;
resampledData = [];
if true%min(pWeights) < maxW
    for s = 1:size(stepVals,2)
            usedIDX = undoneWeights < stepVals(s);
            curData = undoneData(usedIDX,:);
            undoneData(usedIDX,:) = [];
            undoneWeights(usedIDX) = [];
            if not(isempty(curData))      
                if size(curData,1) <= minRes
                    resampledData = [resampledData; curData];
                else
                    n = ceil(size(curData,1)*redRate(s));
                    if n < minRes
                        if meth == 1  
                            %display('resamplePointDataEQ')
                            resampledData = [resampledData; resamplePointDataEQ(curData,minRes,dim)];
                        elseif meth == 2
                            %old version, used for paper
                            resampledData = [resampledData; resamplePointDataByClusters(curData,minRes)];
                            %following is a newer version 
                            %resampledData = [resampledData; resamplePointDataByClustersV2(curData,minRes)];
                        else
                            resampledData = [resampledData; resamplePointData(curData,minRes,dim)];
                        end
                    else
                        if meth == 1 
                            resampledData = [resampledData; resamplePointDataEQ(curData,n,dim)];
                        elseif meth == 2
                            %old version, used for paper
                            resampledData = [resampledData; resamplePointDataByClusters(curData,n)];
                            %following is a newer version 
                            %resampledData = [resampledData; resamplePointDataByClustersV2(curData,minRes)];
                        else
                            resampledData = [resampledData; resamplePointData(curData,n,dim)];
                        end                            
                    end
                end           
            end
    end
end
resampledData = [resampledData; undoneData];
end
%create neightborhood by sectoring
function [ SearchIndxs,NgbrhdMap,stepSize,origin] = createNeighborhoodMap_local(data,step,NeighborhoodLevel)
%Create Search-Indindices and NeighbothoodMap over parameterspace "data" as
%preprocessing to avoid nXn-Search
display('createNeighborhoodMap(...) ... begin preprocessing ...')
s = size(data);
DataMin = min(data);
origin = DataMin;
DataMax = max(data);
nData = s(1);   %number of datasets
dataDim = s(2); %dimension of param-space

if length(step) == 1 || length(step) == dataDim
    stepSize = (DataMax-DataMin)./step+(DataMax-DataMin)./(step)*0.1;
else
    display('createNeighborhoodMap(...) ... Warning, inconvenient "step"-parameter-size!')
    display(' ... User-definded "step" will be ignored, "step"=10 will be used as constant for all dimensions instead.')
    stepSize = (DataMax-DataMin)/10;
end

display('createNeighborhoodMap(...) ... assign search-indices ...')
SearchIndxs = zeros(s);
for d = 1:nData
    for dim = 1:dataDim
           while data(d,dim) >= DataMin(dim) + SearchIndxs(d,dim)*stepSize(dim)
                SearchIndxs(d,dim) = SearchIndxs(d,dim) + 1;
           end
    end
end

uSI = unique(SearchIndxs,'rows'); 

if length(uSI) < 0.75*prod(step)
    display(['createNeighborhoodMap(...) ... Warning, only ',num2str(length(uSI)),' of ',num2str(prod(step)),' possible sectors are occupied!'])
end

display('createNeighborhoodMap(...) ... create 0-neighborhood ...')
s2 = size(uSI);
nIdxs = s2(1);

[def,location] = ismember(SearchIndxs,unique(SearchIndxs,'rows'),'rows');
clear def;

baseMap = cell(nIdxs,2);
NgbrhdMap = cell(nIdxs,3);
for idx = 1:nIdxs
    baseMap(idx,1) = {uSI(idx,:)};
end
for d = 1:nData
    baseMap(location(d),2) = {[(baseMap{location(d),2}) d]};
end

NgbrhdMap = baseMap;
nPts = zeros(size(NgbrhdMap,1),1);
if not(NeighborhoodLevel == 0)
    display(['createNeighborhoodMap(...) ... create ',num2str(NeighborhoodLevel),'-neighborhood ...'])
    for i = 1:nIdxs
        baseIdx = uSI(i,:);
        nn = (NeighborhoodLevel*2+1)^dataDim-1;
        %% new version with given nIndices        
        NIdxs = getNIdxs(baseIdx,NeighborhoodLevel);
        ptIdxs = zeros(1, 10*1024);
        for j = 1:size(NIdxs,1)
            IDX = NIdxs(j,:);
            lc = zeros(size(uSI,1),1);
            for d = 1:dataDim
                lc = lc+(uSI(:,d)==IDX(d));
            end
            loc = find(lc==dataDim,1);            
            if not(isempty(loc)) && loc > 0
                    addIdxs = baseMap{loc,2};
                    nAdd = size(addIdxs,2);
                    if nPts(i)+nAdd > size(ptIdxs,2)
                        ptIdxs = [ptIdxs zeros(1, nn*nAdd)];
                    end
                    ptIdxs((nPts(i)+1):(nPts(i)+nAdd)) = addIdxs;
                    nPts(i) = nPts(i) + nAdd;                
            end           
        end       
        NgbrhdMap{i,2} = [ptIdxs(1:nPts(i)) baseMap{i,2}];           
  
    end
    NgbrhdMap(:,3) = baseMap(:,2);

end
end
function [ SearchIndxs,NgbrhdMap,stepSize,step,origin] = createNeighborhoodMapAdaptive_local(data,step,NeighborhoodLevel)
%Create Search-Indindices and NeighbothoodMap over parameterspace "data" as
%preprocessing to avoid nXn-Search
display('createNeighborhoodMapAdaptive(...) ... begin preprocessing ...')
s = size(data);
DataMin = min(data);
origin = DataMin;
DataMax = max(data);
nData = s(1);   %number of datasets
dataDim = s(2); %dimension of param-space


if length(step) == 1 || length(step) == dataDim
        stepSize = (DataMax-DataMin)./step+(DataMax-DataMin)./(step)*0.1;
else
        display('createNeighborhoodMapAdaptive(...) ... Warning, inconvenient "step"-parameter-size!')
        display(' ... User-definded "step" will be ignored, "step"=10 will be used as constant for all dimensions instead.')
        stepSize = (DataMax-DataMin)/10;
end

display('createNeighborhoodMapAdaptive(...) ... assign search-indices ...')
SearchIndxs = zeros(s);
for d = 1:nData
        for dim = 1:dataDim
               while data(d,dim) >= DataMin(dim) + SearchIndxs(d,dim)*stepSize(dim)
                    SearchIndxs(d,dim) = SearchIndxs(d,dim) + 1;
               end
        end
end

uSI = unique(SearchIndxs,'rows'); 
if true%dataDim == 2

    ratio = ones(dataDim,1);
    for i = 1:dataDim
        ratio(i) = step(1)/step(i);
    end
    nGoal = prod(step);%sx*sy;
    while length(uSI) < 0.75*nGoal
        step(1) = step(1) + ceil(0.1*step(1));
        for i = 2:dataDim
            step(i) = ceil(step(1)/ratio(i));
        end        
        
        stepSize = (DataMax-DataMin)./step+(DataMax-DataMin)./(step)*0.1;
        SearchIndxs = zeros(s);
        for d = 1:nData
                for dim = 1:dataDim
                       while data(d,dim) >= DataMin(dim) + SearchIndxs(d,dim)*stepSize(dim)
                            SearchIndxs(d,dim) = SearchIndxs(d,dim) + 1;
                       end
                end
        end
        display(['createNeighborhoodMapAdaptive(...) ... now ',num2str(length(uSI)),' of ',num2str(nGoal),' wanted sectors are occupied!'])        
        uSI = unique(SearchIndxs,'rows');
    end
end
toc
display('createNeighborhoodMapAdaptive(...) ... create 0-neighborhood ...')
s2 = size(uSI);
nIdxs = s2(1);

[def,location] = ismember(SearchIndxs,unique(SearchIndxs,'rows'),'rows');
clear def;

baseMap = cell(nIdxs,2);
NgbrhdMap = cell(nIdxs,3);
for idx = 1:nIdxs
    baseMap(idx,1) = {uSI(idx,:)};
    
end
for d = 1:nData
    baseMap(location(d),2) = {[(baseMap{location(d),2}) d]};
end
toc
NgbrhdMap(:,1) = baseMap(:,1);
nPts = zeros(size(NgbrhdMap,1),1);
if not(NeighborhoodLevel == 0)
    display(['createNeighborhoodMapAdaptive(...) ... create ',num2str(NeighborhoodLevel),'-neighborhood ...'])
    for i = 1:nIdxs
        baseIdx = uSI(i,:);        
        nn = (NeighborhoodLevel*2+1)^dataDim-1;
        %% new version with given nIndices        
        NIdxs = getNIdxs(baseIdx,NeighborhoodLevel);
        ptIdxs = zeros(1, 10*1024);
        for j = 1:size(NIdxs,1)
            IDX = NIdxs(j,:);
            lc = zeros(size(uSI,1),1);
            for d = 1:dataDim
                lc = lc+(uSI(:,d)==IDX(d));
            end
            loc = find(lc==dataDim,1);            
            if not(isempty(loc)) && loc > 0
                    addIdxs = baseMap{loc,2};
                    nAdd = size(addIdxs,2);
                    if nPts(i)+nAdd > size(ptIdxs,2)
                        ptIdxs = [ptIdxs zeros(1, nn*nAdd)];
                    end
                    ptIdxs((nPts(i)+1):(nPts(i)+nAdd)) = addIdxs;
                    nPts(i) = nPts(i) + nAdd;                
            end           
        end       
        NgbrhdMap{i,2} = [ptIdxs(1:nPts(i)) baseMap{i,2}];            
    end
    NgbrhdMap(:,3) = baseMap(:,2);
end
end
function [SearchIndxs,NgbrhdMap] = createNeighborhoodMapByGivenSektorization_local(data,sektorizeOrigin ,sektorizeStepSize,NeighborhoodLevel)
display('createNeighborhoodMapByGivenSektorization(...) ... begin preprocessing ...')
s = size(data);
nData = s(1);   %number of datasets
dataDim = s(2); %dimension of param-space
DataMin = sektorizeOrigin;
stepSize = sektorizeStepSize;

display('createNeighborhoodMapByGivenSektorization(...) ... assign search-indices ...')
SearchIndxs = zeros(s);
for d = 1:nData
    for dim = 1:dataDim
           while data(d,dim) >= DataMin(dim) + SearchIndxs(d,dim)*stepSize(dim)
                SearchIndxs(d,dim) = SearchIndxs(d,dim) + 1;
           end
    end
end

uSI = unique(SearchIndxs,'rows'); 

display('createNeighborhoodMapByGivenSektorization(...) ... create 0-neighborhood ...')
s2 = size(uSI);
nIdxs = s2(1);

[def,location] = ismember(SearchIndxs,unique(SearchIndxs,'rows'),'rows');
clear def;


baseMap = cell(nIdxs,2);
NgbrhdMap = cell(nIdxs,3);
for idx = 1:nIdxs
    baseMap(idx,1) = {uSI(idx,:)};
end
for d = 1:nData
    baseMap(location(d),2) = {[(baseMap{location(d),2}) d]};    
end
NgbrhdMap = baseMap;
nPts = zeros(size(NgbrhdMap,1),1);
if not(NeighborhoodLevel == 0)
    display(['createNeighborhoodMapByGivenSektorization(...) ... create ',num2str(NeighborhoodLevel),'-neighborhood ...'])
    for i = 1:nIdxs
        baseIdx = uSI(i,:);
        nn = (NeighborhoodLevel*2+1)^dataDim-1;
        %% new version with given nIndices        
        NIdxs = getNIdxs(baseIdx,NeighborhoodLevel);
        ptIdxs = zeros(1, 10*1024);
        for j = 1:size(NIdxs,1)
            IDX = NIdxs(j,:);
            lc = zeros(size(uSI,1),1);
            for d = 1:dataDim
                lc = lc+(uSI(:,d)==IDX(d));
            end
            loc = find(lc==dataDim,1);            
            if not(isempty(loc)) && loc > 0
                    addIdxs = baseMap{loc,2};
                    nAdd = size(addIdxs,2);
                    if nPts(i)+nAdd > size(ptIdxs,2)
                        ptIdxs = [ptIdxs zeros(1, nn*nAdd)];
                    end
                    ptIdxs((nPts(i)+1):(nPts(i)+nAdd)) = addIdxs;
                    nPts(i) = nPts(i) + nAdd;                
            end           
        end
        NgbrhdMap{i,2} = [ptIdxs(1:nPts(i)) baseMap{i,2}];          
    end
    NgbrhdMap(:,3) = baseMap(:,2); 
end

end
%find border sectors
function [ isBorderSektor ] = getBorderSektors_local(sektorIds)
%find all sectors with less than 8 used neightbors (2D)
sS = size(sektorIds);
isBorderSektor = zeros(sS(1),1);
numN = 3^sS(2) - 1;

for i = 1:sS(1)
    curIdx = sektorIds(i,:);
    n = 0;
    for j = 1:sS(1)
        if not(i == j)
            test = max(abs(curIdx-sektorIds(j,:)));
            if test == 1
                n = n+1;
            end
        end
    end
    if n < numN
        isBorderSektor(i) = 1;
    end
end
end
%get indices of neighboring points
function [NIdxs] = getNIdxs(Idx,NLevel,varargin)
    hasLimits = false;
    limits = [];
    idim = size(Idx,2);    
    if nargin > 2
        hasLimits = true;
        linput = varargin{1};
        if size(linput,1) == 1
            limits = repmat(linput,idim,1);
        elseif size(linput,1) == idim
            limits = linput;
        else
            hasLimits = false;
        end
    end
    
    NIdxs = [];
    if NLevel > 0
        NIdxs = zeros((2*NLevel+1)^idim,idim);
        sepIds = zeros(idim,2*NLevel+1);
        
        llow = Idx - NLevel;
        hhigh = Idx + NLevel;
        
        for i = 1:idim
            sepIds(i,:) = [llow(i):hhigh(i)];
            
            tmp = zeros((2*NLevel+1)^i,1);            
            %tmp = [];
            l = (2*NLevel+1)^(i-1);
            for j = 1:2*NLevel+1
                %tmp = [tmp;ones(l,1)*sepIds(i,j)];
                tmp(((j-1)*l)+1:j*l) = ones(l,1)*sepIds(i,j);
            end
            tmp = repmat(tmp,(2*NLevel+1)^(idim-i),1);
            
            %NIdxs = [NIdxs tmp];
            NIdxs(:,i) = tmp;
        end 
         
        NIdxs(ceil(size(NIdxs,1)/2),:) = [];
        
        if hasLimits
            for i = 1:idim
                NIdxs(NIdxs(:,i)<limits(i,1),:) = [];
                NIdxs(NIdxs(:,i)>limits(i,2),:) = [];
            end
        end

    end
end
%resampling a point set
function [choosenDataPoints] = resamplePointData(dataPoints, newPointCount, dim)
%reduce data by construcing regular Points an choose for each regular point
%a corresponding data point!
%input: datapoints : [dim1 dim2 dim3 {params}]
%newPointCount: new number of points per resample dimension
%dim = 1;

s = size(dataPoints);
if s(1) > newPointCount
    dataMax = [max(dataPoints(:,1)) max(dataPoints(:,2)) max(dataPoints(:,3))];
    dataMin = [min(dataPoints(:,1)) min(dataPoints(:,2)) min(dataPoints(:,3))];
    if newPointCount >= 2^dim%dim^2
        ptsPerDim = floor(newPointCount^(1/(dim)));

        step = (dataMax-dataMin)/ptsPerDim;
        numNPts = ptsPerDim^dim;
        RegPt = zeros(numNPts,3);
        %off = (dataMax-dataMin)*0.5;
        %construct regular points
        for t = 1:numNPts
            %compute the nsteps for each index
            if dim > 2
                w = floor((t-1)/ptsPerDim^2);
                %RegPt(t,3) = (0.5*step(3)+w*step(3));
            else
                w = 0;
                %RegPt(t,3) = (off(3));
            end

            if dim > 1
                v = floor(((t-1)-w*ptsPerDim^2)/ptsPerDim);
                %RegPt(t,2) = (0.5*step(2)+v*step(2));
            else
                v = 0;
                %RegPt(t,2) = (off(2));
            end

            u = (t-1) - (w*ptsPerDim^2 + v*ptsPerDim);
            RegPt(t,1) = (0.5*step(1)+u*step(1));
            if dim == 3
                RegPt(t,2) = (0.5*step(2)+v*step(2));
                RegPt(t,3) = (0.5*step(3)+w*step(3));
            elseif dim == 2
                RegPt(t,2) = (0.5*step(2)+v*step(2));

                w = max(u,v);
                RegPt(t,3) = (0.5*step(3)+w*step(3));
            elseif dim == 1
                RegPt(t,2) = (0.5*step(2)+u*step(2));
                RegPt(t,3) = (0.5*step(3)+u*step(3));            
            end

            RegPt(t,:) = dataMin + RegPt(t,:);
        end 
    else
        ptsPerDim = 1;
        numNPts = 1;
        RegPt = mean(dataPoints(:,1:3));
    end
    idxs = zeros(numNPts,1);
    minD = ones(numNPts,1)*max(dataMax-dataMin)*(1/ptsPerDim);

    for t =1:s(1)
        for pn = 1 : numNPts
            d = norm(RegPt(pn,:)-dataPoints(t,1:3));
            if d < minD(pn)
                idxs(pn) = t;
                minD(pn) = d;
            end
        end
    end
    idxs = idxs(not(idxs == 0));
    idxs = unique(idxs);
    choosenDataPoints = dataPoints(idxs,:);
    
    if isempty(choosenDataPoints)
        display(['resamplePointData(...) ....Warning, resampling to zero!'])
    end
else
    choosenDataPoints = dataPoints;
end
end
function [choosenDataPoints] = resamplePointDataEQ(dataPoints, newPointCount, dim)
%reduce data by construcing regular, equdistant Points and choose for each regular point
%a corresponding data point!
%input: datapoints : [dim1 dim2 dim3 {params}]
%newPointCount: new number of points per resample dimension
%dim = 1;

s = size(dataPoints);
if s(1) > newPointCount
    dataMax = [max(dataPoints(:,1)) max(dataPoints(:,2)) max(dataPoints(:,3))];
    dataMin = [min(dataPoints(:,1)) min(dataPoints(:,2)) min(dataPoints(:,3))];
    dataDiff = dataMax-dataMin;
    if newPointCount >= 2^dim

        step = (prod(dataDiff(1:dim))/newPointCount)^(1/dim);
        ptsPerDim = dataDiff(1:dim)/step;
        
        ptsPerDim = floor(ptsPerDim);
        
        if newPointCount==2^dim
            stepDim = 0.5*dataDiff(1:dim);
        else
            stepDim = dataDiff(1:dim)./(ptsPerDim(1:dim)+1);
        end
        
        if dim == 1

            X = [dataMin(1)+0.5*stepDim:stepDim:dataMax(1)-0.5*stepDim];
            RegPt = [X' repmat(mean(dataPoints(:,2:3)),length(X),1)];
        elseif dim == 3
            [X,Y,Z] = meshgrid(dataMin(1)+0.5*stepDim(1):stepDim(1):dataMax(1)-0.5*stepDim(1),dataMin(2)+0.5*stepDim(2):stepDim(2):dataMax(2)-0.5*stepDim(2),dataMin(3)+0.5*stepDim(3):stepDim(3):dataMax(3)-0.5*stepDim(3));
            
            RegPt=[X(:) Y(:) Z(:)];
        else %%dim == 2
            [X,Y] = meshgrid(dataMin(1)+0.5*stepDim(1):stepDim(1):dataMax(1)-0.5*stepDim(1),dataMin(2)+0.5*stepDim(2):stepDim(2):dataMax(2)-0.5*stepDim(2));
            RegPt=[X(:) Y(:) repmat(mean(dataPoints(:,3)),length(X(:)),1)];
        end
        numNPts = size(RegPt,1);
    else
        ptsPerDim = 1;
        numNPts = 1;
        RegPt = mean(dataPoints(:,1:3));
    end
    idxs = zeros(numNPts,1);
    minD = ones(numNPts,1)*step;
    for t =1:s(1)
        for pn = 1 : numNPts
            d = norm(RegPt(pn,:)-dataPoints(t,1:3));
            if d < minD(pn)
                idxs(pn) = t;
                minD(pn) = d;
            end
        end
    end
    idxs = idxs(not(idxs == 0));
    idxs = unique(idxs);
    choosenDataPoints = dataPoints(idxs,:);
    if isempty(choosenDataPoints)
        display(['resamplePointData(...) ....Warning, resampling to zero!'])
    end
else
    choosenDataPoints = dataPoints;
end
end
function [choosenDataPoints] = resamplePointDataByClusters(dataPoints, newPointCount)
% inserted 2014-01-06 by P. Menzel
% kMeans-Clustering for downsampling, resampling points are
% neares-neightbors of cluster centers, dimension-independent
% more irregular result, less artifacts, higher calculation time
% preserves irregular point destributions
    [PointClusterIndex, CusterCenters] = function_kMeans(dataPoints,10,newPointCount,'kMeans++');
    newPts = zeros(size(CusterCenters));
    for c = 1:size(CusterCenters,1)
        center = CusterCenters(c,:);
        resPts = dataPoints(PointClusterIndex == c,:);
        if size(resPts,1)==0
            newPts(c,:) = NaN;
        elseif size(resPts,1)==1
            newPts(c,:) = resPts;
        else
            minD = norm(resPts(1,:)-center);
            id = 1;
            for p=2:size(resPts,1)
                d = norm(resPts(p,:)-center);
                if d < minD
                    minD = d;
                    id = p;
                end
            end
            newPts(c,:) = resPts(id,:);
        end    

    end
    choosenDataPoints = newPts(isfinite(newPts(:,1)),:);
    if isempty(choosenDataPoints)
        display(['resamplePointDataByClusters(...) ....Warning, resampling to zero!'])
    end    
end
%kMeans clustering
function [PointClusterIndex, CusterCenters ] = function_kMeans(data,maxSteps,varargin)
% implmetation for simple kMeans data-clustering
% author: P. Menzel (Geophysik/IFG/CAU Kiel)
% date: 2014-01-06
% version: 1.0 alpha
% http://de.wikipedia.org/wiki/K-Means-Algorithmus 
% based on Kanungo 2002/2004

debug = false;
%% input parameters
initMethode = 'random';
nClusters = 5;

if nargin > 2
    nClusters = varargin{1};
end
if nargin > 3
    initMethode = varargin{2};
end

%% initialize clustering centers
centers = [];
if (strcmp(initMethode, 'regular') && size(data,2)<=3)
    % regular initial cluster centers     
    nClusters = ceil(nClusters^(1/size(data,2)));
    if nClusters < 2
        nClusters = 2;
    end
    
    if size(data,2)==2
        [X,Y] = meshgrid(min(data(:,1)):(max(data(:,1))-min(data(:,1)))/(nClusters-1):max(data(:,1)), ...
                min(data(:,2)):(max(data(:,2))-min(data(:,2)))/(nClusters-1):max(data(:,2)));
        centersReg = [X(:) Y(:)];
    elseif size(data,2)==3
        [X,Y,Z] = meshgrid(min(data(:,1)):(max(data(:,1))-min(data(:,1)))/(nClusters-1):max(data(:,1)), ...
                min(data(:,2)):(max(data(:,2))-min(data(:,2)))/(nClusters-1):max(data(:,2)),...
                min(data(:,3)):(max(data(:,3))-min(data(:,3)))/(nClusters-1):max(data(:,3)));
        centersReg = [X(:) Y(:) Z(:)];        
    else
        centersReg = [min(data):(max(data)-min(data))/(nClusters-1):max(data)]';
    end    
    centers = centersReg;    
elseif (size(nClusters,1) > 1 && size(nClusters,2) == size(data,2))
    % given initial cluster centers     
    centers = nClusters;
elseif strcmp(initMethode, 'kMeans++')
    %http://de.wikipedia.org/wiki/K-Means-Algorithmus#k-Means.2B.2B
    %based on Kanungo 2002/2004
    centersPP = data(randi([1 size(data,1)],1),:);
    for c = 2:nClusters
        distances = zeros(size(data,1),1);        
        for d = 1:size(data,1)
           ld = norm(data(d,:)-centersPP(1,:));
           for c2 = 2:size(centersPP,1)
               dist = norm(data(d,:)-centersPP(c2,:));
               if dist < ld
                   ld = dist;
               end              
           end
           distances(d) = ld;
        end
        centersPP(c,:) = data(find(distances==max(distances),1),:);
    end
    centers = centersPP;    
else
    %random initialization
    loc = zeros(nClusters,1);
    for c = 1:nClusters
        pos = randi([1 size(data,1)],1);
        while not(isempty(find(loc==pos,1)))
            pos = randi([1 size(data,1)],1);
        end
        loc(c) = pos;
    end
    centersRand = data(loc,:);  
    centers = centersRand;
end

%% k-Means
% abort after maxSteps iterations or with changes < abortEps
abortEps = norm(max(data)-min(data))/100; % changes < 1%
it = 1;
do = true;
oldC = centers;
while it <= maxSteps && do
    PCI = zeros(size(data,1),1);
    %assign data to centers
    for d = 1:size(data,1)
       ld = norm2(data(d,:)-centers(1,:));
       cNr = 1;
       for c2 = 2:size(centers,1)
           dist = norm2(data(d,:)-centers(c2,:));
           if dist < ld
               ld = dist;
               cNr = c2;
           end              
       end
       PCI(d) = cNr;        
    end  
    
    %update centers    
    doUpdate = false;
    maxDiff = abortEps;
    for c = 1:size(centers,1)
        if size(data(PCI == c,:),1) == 1
            centers(c,:) = data(PCI == c,:);
        else
            centers(c,:) = mean(data(PCI == c,:)); 
        end
        diff = norm(oldC(c,:)-centers(c,:));
        if diff > abortEps
            doUpdate = true;
            if diff > maxDiff
                maxDiff = diff;
            end
        end       
    end    
    %abort 
    do = doUpdate;
    oldC = centers;
    it = it+1;  
end

PointClusterIndex = PCI;
CusterCenters = centers;
end
%L2-norm for all row-vectors in a matrix
function [sqDist] = norm2(P1)
    %squared euklidean distance
    sqDist = sum(P1.*P1);
end
%estimate and print progress
function [] = printProgress(doneTime,curI,maxI)
%print progress and estimate computation time
%%=========================================================================
%
%   author:       Dipl.-Geoinf. Peter Menzel
%
%   doneTime    = time since computaion begin (toc)
%   curI        = current step
%   maxI        = max step
%%=========================================================================
doTime = (doneTime/curI)*(maxI-curI);
if floor(maxI/100) == 0
    display([' ... Progress: ',num2str(curI),' / ',num2str(maxI),', ready in ',num2str(doTime),'s.'])
elseif mod(curI,floor(maxI/100)) == 0
    display([' ... Progress: ',num2str(ceil((curI/maxI)*100)),'% after ' ,num2str(doneTime),'( /',num2str(doneTime+doTime),')s. Remaining: approx. ',num2str(doTime),'s.'])
end
end
%==========================================================================

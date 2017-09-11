function [rData,effRR] = Other_Resampling_Methods(dataC,dataP,varargin)
%%=========================================================================
%This function contains the algorithms for different resampling approaches,
%which were used to create comparison results in paper GEO-2015-0220: 
% "CIDRe - a parameter constrained irregular resampling method for scattered point data"
%
%
%   author:       Dipl.-Geoinf. Peter Menzel
%   institution:  Christian-Albrechts-Universitaet zu Kiel, Germany, Department for Geosciences
%   date:         2014
%
%   [rData,effRR] = Other_Resampling_Methods(dataC,dataP,dim,steps,methode,resampling_distance)
%   input values:
%       dataC           = data point coordinates
%       dataP           = data point parameter values
%       dim             = 2d(2) or 3D (3)
%       steps           = step to create reguler gridded resampling points
%       methode         = resampling method ('natural','nearest','linear',...
%                           'mundry','selected','minDist','kMeans')
%       resampling_distance
%                       = search distance
%   output values:
%       rData           = resampled dataset
%       effRR           = reduction rate
%==========================================================================

%% input
% geometry dimension
dim = 2;
if nargin > 2
    dim = (varargin{1});
end
% Gridsteps
step = ones(1,dim);
if nargin > 3
    if length(varargin{2}) == dim
        step = varargin{2};
    else
        step = step*varargin{2}(1);
    end
end

% methode
methode = 'natural';
if nargin > 4
    methode = varargin{3};
end

% resampling distance
rDist = 0.1*norm((max(dataC)-min(dataC)));
if nargin > 5
    rDist = varargin{4};
end

%% create Grid
if dim == 1
    gridPts = min(dataC(:,1)):step(1):max(dataC(:,1));
elseif dim == 2 || size(dataC,2) == 2
    [x,y] = meshgrid(min(dataC(:,1)):step(1):max(dataC(:,1)),min(dataC(:,2)):step(2):max(dataC(:,2)));
    gridPts = [x(:) y(:)];
elseif dim == 3
    [x,y,z] = meshgrid(min(dataC(:,1)):step(1):max(dataC(:,1)),min(dataC(:,2)):step(2):max(dataC(:,2)),min(dataC(:,3)):step(3):max(dataC(:,3)));
    gridPts = [x(:) y(:) z(:)];    
end
size(gridPts,1);
size(dataC,1);
RR = size(gridPts,1)/size(dataC,1);

%% create resampled data
rData = zeros(size(gridPts,1),size(dataC,2)+size(dataP,2));
if strcmp('natural',methode) || strcmp('nearest',methode) || strcmp('linear',methode)
    %interpolation auf die Gitterpunkte
    %dataInterpC = cell(1,nDP);
    if dim == 2
        rData(:,1:2) = gridPts;
        if size(dataC,2) == 3
            dataP = [dataC(:,3) dataP];
        end
        for i = 1:size(dataP,2)
            dataInterp = TriScatteredInterp(dataC(:,1),dataC(:,2),dataP(:,i),methode);
            rData(:,2+i) = dataInterp(rData(:,1),rData(:,2));
        end
        
        rData(~isfinite(rData(:,3)),:)=[];
    elseif dim == 3
        rData(:,1:3) = gridPts;
        for i = 1:nDP
            dataInterp = TriScatteredInterp(dataC(:,1),dataC(:,2),dataP(:,i),methode);
            rData(:,3+i) = dataInterp(rData(:,1),rData(:,2));
        end        
        rData(~isfinite(rData(:,4)),:)=[];        
    end
elseif strcmp('mundry',methode)
    % use methode, based on MUNDRY, 1970 interpolation
    rData(:,1:dim) = gridPts;
    for i = 1:size(dataP,2)
        rData(:,dim+i) = mundryV3(dataC(:,1),dataC(:,2),dataP(:,i),gridPts(:,1),gridPts(:,2),rDist,true);
    end
elseif strcmp('selected',methode)
    % selection based on  nearest grid point
    data = [dataC dataP];
    tic
    for g = 1:size(gridPts,1) 
        tmpData = data;
        gP = gridPts(g,:);
        for d = 1:dim
            tmpData = tmpData(tmpData(:,d)>(gP(d)-rDist),:);
            tmpData = tmpData(tmpData(:,d)<(gP(d)+rDist),:);
        end
        
        if numel(tmpData)>0
            minDist = rDist;
            id = 0;
            for i = 1:size(tmpData,1)
                dist = norm(gP(1:dim)-tmpData(i,1:dim));
                if dist < minDist
                    id = i;
                    minDist = dist;
                end
            end
            if id > 0
               rData(g,:) = tmpData(id,:);
            else
               rData(g,:) = nan; 
            end
        else
            rData(g,:) = nan;
        end
        printProgress(toc,g,size(gridPts,1)) 
    end
    rData(~isfinite(rData(:,1)),:)=[];
    rData = unique(rData,'rows');
elseif strcmp('minDist',methode)     
    % used in paper as validation
    rData =zeros(size(dataC,1),size(dataC,2)+size(dataP,2));
    rData(1,:) = [dataC(1,:) dataP(1,:)];
    usedP = 1;
    distmin = rDist;
    distminq=distmin^2;
    tic
    for is=2:size(dataC,1)       
        dist=(dataC(is,1)-rData(1:usedP,1)).^2+(dataC(is,2)-rData(1:usedP,2)).^2;
        if max(dist < distminq) == 0        
            rData(usedP+1,:)=[dataC(is,:) dataP(is,:)];
            usedP = usedP +1;
        end
    end
    rData = rData(1:usedP,:);      
elseif strcmp('kMeans',methode)  
    %[~,newPts]=function_kMeans(dataC,10,nC,'kMeans++');%kMeans++-Initialization
    [~,newPts]=function_kMeansV2(dataC,10,gridPts);%given initial centers
    %sampling based on initial regular cluster centers
    data = [dataC dataP];
    tic
    for g = 1:size(newPts,1) 
        tmpData = data;
        gP = newPts(g,:);
        for d = 1:dim
            tmpData = tmpData(tmpData(:,d)>(gP(d)-rDist),:);
            tmpData = tmpData(tmpData(:,d)<(gP(d)+rDist),:);
        end
        
        if numel(tmpData)>0
            minDist = rDist;
            id = 0;
            for i = 1:size(tmpData,1)
                dist = norm(gP(1:dim)-tmpData(i,1:dim));
                if dist < minDist
                    id = i;
                    minDist = dist;
                end
            end
            if id > 0
               rData(g,:) = tmpData(id,:);
            else
               rData(g,:) = nan; 
            end
        else
            rData(g,:) = nan;
        end
        printProgress(toc,g,size(gridPts,1)) 
    end
    rData(~isfinite(rData(:,1)),:)=[];
    rData = unique(rData,'rows');    
else
    display(['GriddedResampling(...): Warning - "',methode,'"-interpolation is unknown!'])
    rData = [dataC dataP];
end

%% effective resampling rate
effRR = size(rData,1)/size(dataC,1);
display(['GriddedResampling(...): INFO - effective resampling rate: ',num2str((1-effRR)*100),'%'])
end
%% local functions
function ZNEW=mundryV3(X,Y,Z,XNEW,YNEW,minDist,oP)
% reimplementation of FORTRAN-code provided by S. Schmidt, CAU Kiel,
% Department of Geosciences
% original german header:
% C UNTERPROGRAMM ZUR INTERPOLATION BELIEBIG VERTEILTER DATEN (X, Y, Z)
% C AUF BELIEBIG VERTEILTE PUNKTE (XNEW, YNEW) 
% C Methode: MUNDRY, 1970 
% C
% C X,Y         : Koordinaten der originalen Punkte (Vektoren)
% C Z           : Funktionswerte der originalen Punkte (Vektoren)
% C XNEW,YNEW   : Koordinaten der zu interpolierenden Punkte (Vektoren oder 2D Matrizen)
% C ZNEW        : Funktionswerte der zu interpolierenden Punkte (Typ identisch mit XNEW bzw. YNEW)
%
% english header:
% C subroutine for interpolation of scattered data(X, Y, Z)
% C to other scattered locations (XNEW, YNEW) 
% C algorithm based on: MUNDRY, 1970
% C
% C X,Y         : given coordinates (vectores)
% C Z           : given parameter values (Vektoren)
% C XNEW,YNEW   : new locations (vectores oder 2D grids)
% C ZNEW        : interpolated values (type indentically to XNEW bzw. YNEW)
%%
%%=========================================================================
%   change to minDistance mundry, R now has a meaning
%   author:       Dipl.-Geoinf. Peter Menzel (CAU-Kiel/IfG/Abt. Geophysik)
%   date:         2013-01-30
%   version:      0.1 alpha
%
%   update 2014-11-06:  
%                   -minDist===will be handled by global interpolation
%                   -R will be adaptivly increased, when it's to less points
%                       are found
%%========================================================================
%%
n1=size(XNEW,1);
n2=size(XNEW,2);
XNEW=XNEW(:);
YNEW=YNEW(:);
ZNEW=XNEW*0;
N=length(X);

if minDist==0                                                       %P
    gR=1.E+20;                                                      %P
else                                                                %P
    gR = minDist*minDist;                                           %P
end                                                                 %P
if oP
    tic                                                             %Peter
end
for J=1:length(XNEW)
    XX=XNEW(J);
    YY=YNEW(J);
    WMIN=1.E+20;
    W = zeros(N,1);                                                 %P
    for L=1:N
        W(L)=(XX-X(L))^2+(YY-Y(L))^2;
        if(W(L) < WMIN)
            LL=L;
            WMIN=W(L);
        end
    end
    if(WMIN > 1.E-12)
        %R=1.E+20;                                                  %P
        %R = minDist*minDist;                                       %P
        R = gR;                                                     %P
        nN = 0;                                                     %P
        while nN < 6                                                %P
        nN = 0;                                                     %P
        R=2*R;                                                      %P
        S00=0;
        S10=0;
        S01=0;
        S20=0;
        S11=0;
        S02=0;
        Z00=0;
        Z10=0;
        Z01=0;
        for L=1:N
            %if(L ~= LL)%P
            if(L ~= LL && W(L) <= R)                                %P
                %GG=1./W(L)^2;                                      %P
                GG=(R-W(L))/W(L);                                   %P
                S00=S00+GG;
                D1=Z(L);
                U=X(L)-XX;
                V=Y(L)-YY;
                S=GG*U;
                S10=S10+S;
                S11=S11+S*V;
                Z10=Z10+S*D1;
                S20=S20+S*U;
                S=GG*V;
                S01=S01+S;
                S02=S02+S*V;
                Z01=Z01+S*D1;
                Z00=Z00+GG*D1;
                nN=nN+1;                                            %P
            end
        end
        end                                                         %P
        %W1=1/W(LL)^2;                                              %P
        W1 = (R-W(LL)/W(LL));                                       %P
        U=X(LL)-XX;
        V=Y(LL)-YY;
        Z1=Z(LL);
        if(abs(U) > abs(V))
            A11=U*S00-S10;
            A21=U*S10-S20;
            A31=U*S01-S11;
            GG=1/(W1*U);
            A12=1.+GG*S10;
            A22=U+GG*S20;
            A32=V+GG*S11;
            GG=V/U;
            A13=S01-GG*S10;
            A23=S11-GG*S20;
            A33=S02-GG*S11;
            R1=U*Z00-Z1*S10;
            R2=U*Z10-Z1*S20;
            R3=U*Z01-Z1*S11;
        else
            A11=V*S00-S01;
            A21=V*S10-S11;
            A31=V*S01-S02;
            GG=1./(W1*V);
            A12=1.+GG*S01;
            A22=U+GG*S11;
            A32=V+GG*S02;
            GG=U/V;
            A13=-S10+GG*S01;
            A23=-S20+GG*S11;
            A33=-S11+GG*S02;
            R1=V*Z00-Z1*S01;
            R2=V*Z10-Z1*S11;
            R3=V*Z01-Z1*S02;
        end
        
        D1=A22*A33-A32*A23;
        D2=A12*A33-A32*A13;
        D3=A12*A23-A22*A13;
        DET=A11*D1-A21*D2+A31*D3;
        if DET~=0 
            ZNEW(J)=(R1*D1-R2*D2+R3*D3)/DET;
        end
    else
        ZNEW(J)=Z(LL,1);
    end
    if oP
        printProgress(toc,J,length(XNEW))                             %Peter
    end
end

ZNEW=reshape(ZNEW,n1,n2);
end
function [ PointClusterIndex, CusterCenters ] = function_kMeansV2(data,maxSteps,varargin)
% implmetation for simple kMeans data-clustering
% author: P. Menzel (Geophysik/IFG/CAU Kiel)
% date: 2014-11-04
% version: 1.0 alpha
% based on Knungo et. al. 2002 and 2004
% optimized Version

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
display(['... kMeans: Initialize -',initMethode,'-'])
if (strcmp(initMethode, 'regular') && size(data,2)<=3)
    % initial regularly distributed cluster centers       
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
    % given cluster centers
    centers = nClusters;
elseif strcmp(initMethode, 'kMeans++')
    %http://de.wikipedia.org/wiki/K-Means-Algorithmus#k-Means.2B.2B
    %Kanungo et. al. 2002/2004
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
% abbort after maxSteps iterationens or no further changes in the cluster
% centers
display(['... kMeans: Clustering '])
abortEps = norm(max(data)-min(data))/100; % update difference < 1%
it = 1;
do = true;
oldC = centers;
while it <= maxSteps && do
    PCI = zeros(size(data,1),1);
    %assignement points to centers    
    display(['....... assign centers'])
    dMat = getNormVec(data,centers(1,:),0);
    PCI(:)=1;
    for c = 2:size(centers,1)
        tmp = getNormVec(data,centers(c,:),0);
        loc = tmp<dMat;
        PCI(loc) = c;
        dMat(loc) = tmp(loc);
    end    
    display(['....... update centers'])
    %update der centers    
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
    %abort stuff  
    do = doUpdate;
    oldC = centers;
    display(['... done Iteration ',num2str(it),'/',num2str(maxSteps),' - maxDiff = ',num2str(maxDiff)])
    it = it+1;  
end
%% return 
PointClusterIndex = PCI;
CusterCenters = centers;
end
function [nVec] = getNormVec(m,v,squared)
%Matrix m (nxv): consists of a number of row-vectors,
%v: single row vector (1xv), 
%nVec: vector (nx1) with l2-Norm of rows with |m(i,:)-v|
s = size(m);
nVec = zeros(s(1),1);
for d = 1:s(2)
    dVec = m(:,d)-v(d);
    nVec = nVec+dVec.*dVec;
end
if squared
    nVec=sqrt(nVec);
end
end
function [] = printProgress(doneTime,curI,maxI,varargin)
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
    display([' ... progress: ',num2str(curI),' / ',num2str(maxI),', ready in ',num2str(doTime),'s.'])
elseif nargin == 3 && mod(curI,floor(maxI/100)) == 0
    display([' ... progress: ',num2str(ceil((curI/maxI)*100)),'% after ' ,num2str(doneTime),'( /',num2str(doneTime+doTime),')s. Remaining: approx. ',num2str(doTime), 's.'])
elseif nargin > 3 && mod(ceil((curI/maxI)*100),varargin{1}) == 0
    display([' ... progress: ',num2str(ceil((curI/maxI)*100)),'% after ' ,num2str(doneTime),'( /',num2str(doneTime+doTime),')s. Remaining: approx. ',num2str(doTime), 's.'])
end
end
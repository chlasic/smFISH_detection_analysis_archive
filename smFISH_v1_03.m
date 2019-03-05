%%% This code recognizes and records coordinates of nuclei and transcripts in each germline image in a folder. 
%% Clear workspace and create folders to save outputs.
% clc; clearvars; 
tic;

%%% Settings for image processing and analysis
%%% Create output folder and subfolders

%%%%%%------------- User Input (UI) system --------------------
% fpath = uigetdir('D:\smFISH images\');
% temp = regexp(fpath, '\', 'split');
% IMAGEfolderName = temp{end}; clear temp;
%-----------------------------------------------------------

%%%%%%------------ Manual input -------------------------------
%%% -------  for Windows
MASTER_Image_path = 'C:\CHL_N2'; % <<<<<-- Change the path!!!  Master image folder (one folder above the folder containing the image file to process).
IMAGEfolderName   =  '';           % <<<<<-- Change the path!!!  Input folder containing the image file (folder below the Master image folder).

MASTER_Output_path = 'D:\smFISH analyses\'; % <<<<<-- Change the path!!!  Master Output folder
OUTPUTdir1=  'CHL_new_N2';
% OUTPUTdir1=  strcat(IMAGEfolderName,'_result');

% Create output folder and subfolders
MASTER_Image_path = strcat(MASTER_Image_path, IMAGEfolderName, '\');
cd(MASTER_Output_path); mkdir(OUTPUTdir1); cd(OUTPUTdir1);
savepath = strcat(MASTER_Output_path,OUTPUTdir1,'\');
%%%-----------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Threshold for image processing (-0.5 to inf; normally 0 - 1). Change parameter if not detect all spots
% Thrshold for Exon/Intron:
% N2: .5/.5     lag-2(q411)/+: .5/.7     glp-1(q46)/+: .5/.7     lag-1(q385)/+: .3/.7   JK5008: .65/.45 .
% lag-3(ok387) het: .5/.4
% q224: 
% Erika's mex-3: .5/.5
thresForExon =         .2     ; 
thresForIntron =      1.1     ;
thresForNuc =    1      ;

nthfile =         1      ;       

% Designate the channel for RNA.
% Put   '1' for intron,     '2' for exon,      '3' for else.
ChOrder = [     2       1        ];


radius =     2.5     ;   % set the Radius of ROI as desired (um)
vLim =      3        ;   % maxium distance for Voronoi cell.
nrange = [   1.2     2.8   ];    % define range of nuclear radius in um.
sensi =     0.96   ;     % sensitivity for nuclear circle detection (higher: more flex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% %%% --------  for MAC --------------------------------------                                                                                         
% MASTER_Image_path = strcat('/Users/hjshin/Documents/MATLAB/smFISH/'); % <<<<<-- Change the path!!!  Master image folder (one folder above the folder containing the image file to process).
% IMAGEfolderName   =  'sel-12';                   % <<<<<-- Change the path!!!  Input folder containing the image file (folder below the Master image folder).
% 
% MASTER_Output_path=  '/Users/hjshin/Documents/MATLAB/smFISH/';              % <<<<<-- Change the path!!!  Output folder
% OUTPUTdir1=  'sel-12 condition 1';
% % OUTPUTdir1=  IMAGEfolderName;
% 
% % Create output folder and subfolders
% MASTER_Image_path = strcat(MASTER_Image_path, IMAGEfolderName, '/');
% cd(MASTER_Output_path); mkdir(OUTPUTdir1); cd(OUTPUTdir1);
% savepath = strcat(MASTER_Output_path,OUTPUTdir1,'/');
% %%%-----------------------------------------------------------  


%%% Open matlabpool for parallel computing
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

fprintf('\n\t\tThe radius of ROI is %2.1f um.\n\n', radius); 


%% Read lif files 
%%% lif file format: separate diff condition by ' ' (space). Order or # of
% statements does not matter. e.g. 040914 rrf-1 emptyRNAi_na lag-3_na.lif.
% liflist = row1: lif file name w/o '.lif'; row2: name with sorted; conditions; row3: strain numbering.
[liflist] = dir(strcat(MASTER_Image_path,'*.lif'));
liflist = liflist(~[liflist.isdir]);
liflist = {liflist.name}';
nbf = numel(liflist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Choose n-th lif file to process
nbf = 1;
liflisttemp = liflist(          nthfile           );
liflist=liflisttemp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% gather lif files with the same labels (experimental conditions).
liflist = strrep(liflist,'.lif','');
lifSplit = regexp(liflist, ' ', 'split');

parfor i = 1:nbf
    lifSplit{i} = lifSplit{i}(cellfun('isempty', regexp(lifSplit{i},'^-?\d+$')));
    lifSplit{i} = sort(lifSplit{i}); 
    liflist{i,2} = strjoin(lifSplit{i},'_');
end

[listSort,indeX] = sortrows(liflist, 2);
listSort{1,3} = 1; liflist{indeX(1),3} = 1;

count = 1;
st_names = cell(nbf,1);
st_names{1} = listSort{1,2};
for i = 1:nbf-1
    if strcmp(listSort{i+1,2},listSort{i,2}) == 1
        listSort{i+1,3} = listSort{i,3};
        liflist{indeX(i+1),3} = liflist{indeX(i),3};
    else
        count = count+1;
        listSort{i+1,3} = count;
        liflist{indeX(i+1),3} = count;
        st_names{count,1} = listSort{i+1,2};
    end
end

if nbf ~= 1
    st_names(count+1,:) = [];
end

clear count
st = length(st_names);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stname_forNum = cell(1,listSort{end,end}); 
af = cell(999,10,nbf);


%% Process multiple lif files sequentially
for h=1:nbf    % h is a counter for lif files (with diff conditions) in the folder.
%     clearvars -except af pix sensi nrange st_names indeX liflist h listSort pDir1...
%         stname_forNum CIpercent nbf radius label1 label2 label3 label4 fpath savepath...
%         thresForExon thresForIntron thresForNuc ChOrder vLim kkk kno nthfile ifo;
        
%%% import lif files in the same conditions (e.g. strain, temp)
lifdat = cell(1,4);
fToRead = strcat(MASTER_Image_path, listSort{h,1}, '.lif');

% read in total # of image stacks
lifRead = bfGetReader(fToRead);
omeMeta = lifRead.getMetadataStore();
lic = omeMeta.getImageCount();


    %% Process and analyze individual images one by one.
%     kk = 4:5;
    for f =  1:lic     % open f-th image-stack in the lif file.
        
%         if f == 19
%             continue
%         end
        
%         if f == 4
%             thresForNuc =    .8      ;
%         elseif f == 5
%             thresForNuc =    1      ;
%         else
%             thresForNuc =    .5      ;
%         end

%         if f == 9 || f == 13 || f == 19
%             thresForExon =         .5     ; 
%         end
        
        lifdat = bfopenOne(fToRead,f);

        %%% read metadata
        meta = lifdat{1, 4};
        Nch = meta.getPixelsSizeC(f-1).getValue();
        voxelSizeXdefaultValue = meta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
        voxelSizeXdefaultUnit = meta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
        voxelSizeX = meta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); 
        voxelSizeY = meta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); 
        voxelSizeZ = meta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER); 
        
        iminfo = [
        Nch % # channel
        meta.getPixelsSizeX(f-1).getValue(); % image width, pixels
        meta.getPixelsSizeY(f-1).getValue(); % image height, pixels
        meta.getPixelsSizeZ(f-1).getValue(); % number of Z slices
        voxelSizeX.doubleValue();  % in µm
        voxelSizeY.doubleValue();  % in µm
        voxelSizeZ.doubleValue()   % in µm
        ];
                                      

        %%% in case above method does not work for reading # of z-planes-----
%         temp = meta.getPixelsPhysicalSizeZ(f-1).getValue(); % in um, works on ver.R2014b and +.
%         if isnan(temp) == 1
%             str2double(lifdat{1,2}.get('ATLConfocalSettingDefinition|Quantity|Value 0')) * 1000000;
%         end
%         iminfo = [iminfo; temp];
        %---------------------------------------------------------------------

        %%% Import image information from metadata
        % find Nuc, phase or trx channel by color (channel label does not work).  
        % Chinfo: |col 1: color info in the order of acquisition|  
        Chinfo = zeros(Nch,2);
        for i = 1:Nch
            Chinfo(i,1) = meta.getChannelColor(f-1,i-1).getValue();
        end

        % Nuc | phase | Yellow (often 561/594) | Red (often 633/594) | green (often LMN-1 or other Ab with Alexa488)
        % 1: DAPI(blue or cyan), 2: Phase, 3: Yellow, 4: Red, 5: Magenta, 6: Green
        Chinfo(Chinfo(:,1) == 16777215 | Chinfo(:,1) == 65535,2) = 1; % DAPI, blue or cyan
        Chinfo(Chinfo(:,1) == -1,2) = 2; % gray, NOT phase channel
        Chinfo(Chinfo(:,1) == -65281,2) = 3; % yellow, 1st RNA
        Chinfo(Chinfo(:,1) == -16711681,2) = 4; % magenta, 2nd RNA?
        Chinfo(Chinfo(:,1) == -16776961,2) = 5; % red, 3rd RNA     
        Chinfo(Chinfo(:,1) == 16711935,2) = 6; % green, Ab staining? (e.g. Alexa-488)
        if find (Chinfo(:,1) == -1, 1, 'last') == Nch
            Chinfo(Nch,2) = 9; % phase
        end

        loc = cell2mat(strfind(lifdat{1,1}(:,2), 'C='))+2;
        chameta = lifdat{1,1}(:,2);
        leng = length(chameta(:,1));
        chainfo = zeros(leng,1);
        parfor i = 1:leng
            temp = char(chameta(i,1));
            chainfo(i) = str2double(temp(loc(i)));
        end



        

        %% Detect a germline outline using ATS channel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Seperate planes
        % The last number on line below determines which color channel to use.
        PhaPla = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == 3), :);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pix = round( 0.22 / mean(iminfo(5:6,1))); 
        if pix == 0
            pix = 1;
        end

        zproj = cat(3, PhaPla{:,1});
        mip = max(zproj, [], 3);
        mipBW = imbinarize(mip, graythresh(mip)/3);


        %%% "cyl_axis" contains: (unit: pixels)
        %%% (1) x-coor | (2) y-coor | (3) z-coor | (4) radius of the cell on that point.| 
        [cyl_axis, germlineBound,cmask] = DefOutlineATS(mip, pix, iminfo);

        %%% align germline orientation (distal end to the left)
        if mean(cyl_axis(1:10,4)) > mean(cyl_axis(end-9:end,4))
            cyl_axis(:,1) = iminfo(2)-cyl_axis(:,1)+1;
            cmask = flip(cmask, 2);
            lifdat{1,1}(:,1) = cellfun(@(x) flip(x,2), lifdat{1,1}(:,1), 'uni', 0);
        end


        %% Detect Nuclei

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % seperate Nuc channel to detect nuclei.
        NucPla = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == 1), :);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Detect nuclei and remove one nucleus if two are overlapped.
        % NDNuc saves not overlapped Nuclei that are in DetNuc (detected Nuc).
        tnran = round(nrange ./ iminfo(5:6,1)'); 
        tol = round( 1.2 / mean(iminfo(5:6,1))); %%%%% tolerance for overlapped circles (um).

        % Detect nuclei from image slices and make nuclear sphere. 
         
        
        [Nuc, thrNuc, aftNuc, tickn] = DetectNucleus(NucPla, iminfo, pix, tnran, sensi, tol, f, lic, thresForNuc);

        % Detect nuclei using bwconncomp and cross-validate the result with nuclei detected by hough transf.
        % 'nucs' : | x-coor | y | z | radius | ave. circularity |
        %           | std. circularity | total DAPI sig. area (# pixel) | 
        nucs = DetectNucBlobs(thrNuc, aftNuc, Nuc, iminfo, tnran);

        % Remove nuclei that are outside the gonad
        nucs(:,8) = 1;
        for i = 1:size(nucs,1)
            r1 = round(nucs(i,2))-3;
            r2 = round(nucs(i,2))+3;
            r3 = round(nucs(i,1))-3;
            r4 = round(nucs(i,1))+3;
            if r1 < 1
                r1 = 1;
            end
            if r2 > iminfo(3)
                r2 = iminfo(3);
            end
            if r3 < 1
                r3 = 1;
            end
            if r4 > iminfo(2)
                r4 = iminfo(2);
            end
            if sum(sum(mipBW(r1:r2, r3:r4))) == 0
                nucs(i,8) = 0;
            end
        end
        nucs(nucs(:,8) == 0,:) = [];
        nucs(:,8) = [];

        % Visualize of nuclei
        %----------- max projection  ------------------
        zproj = cat(3, NucPla{:,1});
        mip = max(zproj, [], 3);
%         figure, imshow(mip);


        %% Detect RNA spots
        % Process different RNA (channels) separately.
        rch = [3 4 5 6];
        chname = {'yel', 'mag', 'red', 'grn'};
        sele = ismember(rch, Chinfo(:,2));
        chnn = rch(sele); 
        
        chn = zeros(1,6);
        loc = 1;
        tempc = Chinfo(:,2);
        for i = 1:length(tempc)
            if ismember(tempc(i),rch)
                chn(loc) = tempc(i);
                loc = loc + 1;
            end
        end
        chn(loc:end) = [];
        
        nch = length(chn);
        chname = chname(sele);
        prna = cell(nch,1);
        nuc = cell(nch,1);
        
%         nch = length(ChOrder);
        
        
        for i = 1:nch
            rnaPla.(chname{i}) = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == chn(i)), :);
            if i == 1
                k = 10; % 1st RNA (channel) image is saved in 10th column
            elseif i == 2
                k = 11; % 2nd RNA (channel) image is saved in 11th column
            else
                k = 12; % 3rd RNA (channel) image is saved in 12th column
            end
            zproj = cat(3, rnaPla.(chname{i}){:,1});
            zproj = max(zproj, [], 3);
            af{f,k} = zproj;
        end

        %%% process channels seperately for detecting RNA
        % 'blobrna': | rna ID | x-size | y-size | z-size | 5th: total # pixel | 
        %            | total intensity | x centroid | y centroid | z centroid.

        % blrna: |Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
        %        | 6th: matching nuc #(row # in 'nuc') 

        % 'rnas' 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
        %               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
        %               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
        %               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|

        % 'nuc for ATS' : | x-coor | y | z | radius | ave. circularity |
        %         | 6th: std. circularity | total DAPI sig. | of all nuclei | # trx sites | RNA ID .

        % 'nuc for mRNA' : | 1:3 xyz-coor | 4: radius | 5: ave. circularity |
        %           | 6: std. circularity | 7: total DAPI sig. area (# pixel) | 8: nuc ID.
        %           | 9: Voronoi mRNA count | 10: Voronoi count w/ limit (vLim, 3um) | 11: # mRNA in ROI w/ radius (2.5um)  |
        %  (aft_analysis->)       | 12th col: dist betw. nuc & cell outline | 13: 0 = cells on cortex, 1 = cells in rachis |

        pixr = round( 0.1 / mean(iminfo(5:6,1)));

        parfor g = 1:nch
            rPla = rnaPla.(chname{g})(:,1);
            %%%%%%%%%%%%%% Detect RNA spots %%%%%%%%%%%%%%%%%%%%%%%%%
            %-------------- process INTRON -----------------------------
            if ChOrder(g) == 1    % 1 is for processing intron images.
                % mom1 or mom2 is the mean background signal intensity from all z-planes in a gonad.
                [blobrna, blrna, mom1(g)] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, lic, g, thresForIntron);

                if isempty(blobrna) == 0
                    % determine which nuclei RNA spots belong. 
                    [prnas, nuct] = MatchrnaNuc(thrNuc, blobrna, nucs, iminfo, pixr, 2);

                    % Put results of RNA1 and RNA2 in separate struct.
                    prna{g} = prnas;
                    nuc{g} = nuct;
                else
                    prna{g} = [];
                    nuc{g} = [];
                end

            % ------------- process EXON image ----------------------------    
            elseif ChOrder(g) == 2   % 2 is for processing exon images.
                [DETrna, derna, mom2(g)] = DetectRNAExon(rPla, pixr, iminfo, f, lic, g, thresForExon);

                if isempty(DETrna) == 0

                    % Put results of RNA1 and RNA2 in separate struct.
                    DETrna(:,8) = 1:length(DETrna(:,1)); % mRNA ID.
                    [mrnas, rnas_atsEx, nucm, nuc_atsEx] = MatchmrnaNuc(DETrna, thrNuc, nucs, iminfo, radius, vLim, pixr);

                    prna{g} = mrnas; % save info for detected mRNA
                    prna{g}(:,4) = 1;
                    nuc{g} = nucm;
                    
                    %%% record detected ATS from exon channel
                    ExATSrna{g} = rnas_atsEx;
                    ExATSnuc{g} = nuc_atsEx;
                else
                    prna{g} = [];
                    nuc{g} = [];
                end
            elseif ChOrder(chn(g)-2) == 3   % 3 is for processing etc.

            end

        end 



        %% cross-validate trx sites by Intron probe and by Exon probe.
        % Also remove mRNA spots that are determined as ATS.
        
        % 'DexRext' 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
        %               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
        %               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
        %               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|

        if ismember(2, ChOrder)
            if ~isempty(prna{ChOrder == 2})
                mRNAIntM = mean(prna{ChOrder == 2}(:,6)); % mean intensity of all mRNA spots.
            else
                mRNAIntM = 0;
            end
        else
            mRNAIntM = 0;
        end

        if ismember(1, ChOrder)
            if ~isempty(prna{ChOrder == 1})
                pRNAIntM = mean(prna{ChOrder == 1}(:,6)); 
            else
                pRNAIntM = 0;
            end
        else
            pRNAIntM = 0;
        end

        if ismember(2, ChOrder) && ismember(1, ChOrder)
            prna{ChOrder == 1}(:,8:9) = 0;
            prna{ChOrder == 2}(:,13) = 0;
            distATS =       1.5            ;   % acceptible distance (um) of ATSs on different channels.

            countp = 0;
            countm = 0;

            % remove detected ATS spots from the list of mRNAs
            if ~isempty(prna{ChOrder == 1})
                for i = 1:length(prna{ChOrder == 1}(:,1))
                    zpln = prna{ChOrder == 1}(i,7);
                    DetRext = prna{ChOrder == 2}(prna{ChOrder == 2}(:,7) > zpln- 6 & prna{ChOrder == 2}(:,7) < zpln+ 6,:);
                    distMat = zeros(length(DetRext(:,1)),5);
                    distMat(:,1) = DetRext(:,1) - prna{ChOrder == 1}(i,1);   
                    distMat(:,2) = DetRext(:,2) - prna{ChOrder == 1}(i,2);  
                    distMat(:,3) = DetRext(:,6);
    %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder == 1}(i,3);   
                    distMat(:,4) = sqrt(distMat(:,1).^2 + distMat(:,2).^2); % + distMat(:,3).^2);
                    distMat(:,5) = DetRext(:,8);
                    distMats = distMat(distMat(:,4) <= distATS/iminfo(6), :);
                    if isempty(distMats)
                        prna{ChOrder == 1}(i,9) = 999;
                        continue
                    end
                    distMats = sortrows(distMats,-3);

                    if distMats(1,4) <= distATS/iminfo(6)
                        if distMats(1,3) < mRNAIntM/2 && prna{ChOrder == 1}(i,6) < pRNAIntM/2 && prna{ChOrder == 1}(i,1) > min(prna{ChOrder == 1}(:,1)) + 25/iminfo(6)
                            prna{ChOrder == 1}(i,9) = 999;  % mark intron spot unlikely ATS.
                            countp = countp + 1;
                        else
                            prna{ChOrder == 2}(prna{ChOrder == 2}(:,8) == distMats(1,5),13) = 999;  % mark ATS in exon channel
                            prna{ChOrder == 1}(i,8) = distMats(1,3);
                            countm = countm + 1;
                        end
                    else
                        prna{ChOrder == 1}(i,9) = 999;
                    end
                end
            end

            % remove intron ATS that are not seen at exon channel
    %         NoOvlapATS = prna{ChOrder == 1}(prna{ChOrder == 1}(:,8) == 999,:);
            if countp > 0
%                 prna{ChOrder == 1}(prna{ChOrder == 1}(:,9) == 999,:) = [];
            end
            prna{ChOrder == 1}(:,9) = [];

            if countm > 0
%                 prna{ChOrder == 2}(prna{ChOrder == 2}(:,13) == 999,:) = [];  
            end
            prna{ChOrder == 2}(:,13) = [];
        end


        
        %%%%%     Record ATS intensity    %%%%%
        %%% 'nuc' for ATS
        %       |1-3: XYZ-coordinates |4: radius (pixels)| 5: ave. circularity| 6: std circularity |
        %       | 7: DAPI intensity |8: # ATS per nuc | 9: summed ATS intensity from INTRON channel |
        %       | 10: summed ATS int. from EXON channel
        if ~isempty(nuc{ChOrder == 1})
            if length(ChOrder) > 1
                nNcount = find( nuc{ChOrder == 1}(:,8) > 0 );
                for i = 1:length(nNcount)
                    ints = prna{ChOrder == 1}(prna{ChOrder == 1}(:,4) == nNcount(i),6);
                    nuc{ChOrder == 1}(nNcount(i),9) = sum(ints);
                    ints2 = prna{ChOrder == 1}(prna{ChOrder == 1}(:,4) == nNcount(i),8);
                    nuc{ChOrder == 1}(nNcount(i),10) = sum(ints2);
                end
            end
        end
        
        

        %% Make summary data
        % 'af': | 1st: nuc1 info | 2nd: rna1 info (yellow) | 3rd: nuc2 info | 4th: rna2 info (red) | .
        %           | 5th: cyl_axis info | 6th: cell outline mask | 7th: z-prj nuclear channel |.
        %           | 8th: mean bg intensity from 1st channel | 9th: mean bg. noise from 2nd channel|
        %           | 10th: z-projected RNA image (1st channel) | 11th: z-prj. 2nd RNA channel |
        %  (aft_analysis->)          | 12th col: distance distal tip & first chain nuc | 12: # nuc in rachis |
        for i = 1:nch
            if ~isempty(prna{i})
                temp = prna{i}(prna{i}(:,4) ~= 0,:);
                af(f,i*2-1, nbf) = nuc(i);
                af(f,i*2, nbf) = {temp};
            else
                af(f,i*2-1, nbf) = {[]};
                af(f,i*2, nbf) = {[]};
            end
        end

        af{f,5, nbf} = cyl_axis;
        af{f,6, nbf} = cmask;
        af{f,7, nbf} = mip;
        if exist('mom1', 'var')
            af{f,8, nbf} = mom1(mom1~=0);
        end
        if exist('mom2', 'var')
            af{f,9, nbf} = mom2(mom2~=0);
        end

        %%% record information for nuc (12 col) & detected ATS (13th col) from exon channel in 'af'
        if exist('ExATSrna', 'var')
            if ismember(2, ChOrder)
                for i = 1:size(ExATSrna,1)
                    if ~isempty(ExATSrna{i})
                        temp = ExATSrna{i}(ExATSrna{i}(:,4) ~= 0,:);
                        af(f,i*2-1 +11, nbf) = ExATSnuc(i);
                        af(f,i*2 +11, nbf) = {temp};
                    else
                        af(f,i*2-1 +11, nbf) = {[]};
                        af(f,i*2 +11, nbf) = {[]};
                    end
                end
            end
        end
        
        af{f,16,nbf} = iminfo;
        
        cd(savepath);
        save(strcat('smFISHworkspace_middle_',num2str(f)), '-regexp', '^(?!(lifdat|lifRead|meta|omeMeta)$).');


    end

    if nbf == 1
        af(lic+1:end,:) = [];
    else
        af(lic+1:end,:,nbf) = [];
    end
  
    
end

cd(savepath);
save('smFISHworkspace_final', '-regexp', '^(?!(lifdat|lifRead|meta|omeMeta)$).');

eTime = toc;
fprintf('\n\t\tElapsed time is %2.1f minutes.\n\n', eTime/60); 
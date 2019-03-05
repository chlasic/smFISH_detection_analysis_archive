function [rnas, rnas_atsEx, nuc, nuc_atsEx] = MatchmrnaNuc(DETrna, thrNuc, nucs, iminfo, radius, vLim, pixr)
% This function matches rnas spots with nuclei (only mRNA spots)

%%% determine which nuc the RNA spots belong. 
% 'rnas' 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
%               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
%               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
%               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|

rnas = DETrna;
rnas_atsEx = rnas;
rnas(:,9:12) = 0;
nuc = nucs;
nuc_atsEx = nuc;
nuc(:,8) = 1:size(nuc,1);
nuc(:,9:11) = 0;
nuc_atsEx(:,8:9) = 0;
rnasNucList = cell(size(rnas,1),1);
sen = strel('disk', pixr+1); 
meanIntensity_mRNA = mean(rnas(:,6));  

% 'nuc' : | 1:3 xyz-coor | 4: radius | 5: ave. circularity |
%           | 6: std. circularity | 7: total DAPI sig. area (# pixel) | 8: nuc ID.
%           | 9: Voronoi mRNA count | 10: Voronoi count w/ limit (vLim, 3um) | 11: # mRNA in ROI w/ radius (2.5um)  |


for i = 1:size(rnas,1)
    cutset = nuc(nuc(:,3) >= rnas(i,3) - radius/iminfo(6) & nuc(:,3) <= rnas(i,3) + radius/iminfo(6),:);
    distemp = zeros(size(cutset,1),2);
    distemp(:,1) = sqrt( (cutset(:,1)-rnas(i,1)).^2 + (cutset(:,2)-rnas(i,2)).^2 + (cutset(:,3)-rnas(i,3)).^2 );
    distemp(:,2) = cutset(:,8);
    distemp(:,3:6) = cutset(:,1:4);
    
    % 'distemp': | distNR (1st col) | nuc ID (2 col) | 
    %            | x-z nuc coor (3-5 col) | nuc radius (6 col) in pixels |
    % select 3 shortest length from nuc center.
    distemp = sortrows(distemp, 1);
    
    if ~isempty(distemp)
        %%% Test if detected exonal ATS spot is on DAPI (overlap).
        if ismember(round(rnas(i,7)), [1 2 3 4]) == 1
            tb = cat(3,thrNuc{1:round(rnas(i,7))+4,1});
        elseif ismember(round(rnas(i,7)), [iminfo(4)-3 iminfo(4)-2 iminfo(4)-1 iminfo(4)]) == 1
            tb = cat(3,thrNuc{round(rnas(i,7))-4:end,1});
        else
            tb = cat(3,thrNuc{round(rnas(i,7))-4:round(rnas(i,7))+4,1});
        end
        tb = max(tb,[],3);
        tx = imdilate(tb,sen);

        tstart = round(rnas(i,2) - pixr*2-1):round(rnas(i,2) + pixr*2+1);
        tend = round(rnas(i,1) - pixr*2-1):round(rnas(i,1) + pixr*2+1);

        tstart(tstart < 1) = [];
        tstart(tstart > iminfo(3)) = [];
        tend(tend < 1) = [];
        tend(tend > iminfo(2)) = [];

        txpart = tx(tstart, tend);
        
        dapiOL = 0;
        if sum(sum(txpart)) >= pixr*           2     %
            dapiOL = 1;
        end
        
        % record exonal ATS seperately.
        if distemp(1,1) < distemp (1,6) && dapiOL == 1
            if rnas_atsEx(i,6) > meanIntensity_mRNA
                rnas_atsEx(i,4) = distemp(1,2);
                nuc_atsEx(distemp(1,2),8) = nuc_atsEx(distemp(1,2),8) + rnas_atsEx(i,5);
                if nuc_atsEx(distemp(1,2),8) > 4
                    nuc_atsEx(distemp(1,2),8) = 4;
                end
                nuc_atsEx(distemp(1,2),9) = nuc_atsEx(distemp(1,2),9) + rnas_atsEx(i,6);
            end
        else
            %%% count mRNAs per cell (ID-ed by nucleus) with 3 diff ways.
            rnas(i,4) = 1;
            % simple Voronoi count
            nuc(distemp(1,2),9) = nuc(distemp(1,2),9) + 1;
            rnas(i,9) = 1;
            rnas(i,12) = distemp(1,1); % unit: pixel

            % Voronoi w/ distance limit
            if distemp(1,1) <= vLim/iminfo(6)
                nuc(distemp(1,2),10) = nuc(distemp(1,2),10) + 1;
                rnas(i,10) = 1;
            end

            % count # mRNA in ROI (rad: Radius; 2.5 um)
            valDist = distemp(distemp(:,1) < radius/iminfo(6),:);
            rnas(i,11) = length(valDist(:,1));
            rnasNucList{i} = [rnasNucList{i} valDist(:,2)'];

            if ~isempty(valDist)
                nuc(valDist(:,2),11) = nuc(valDist(:,2),11) + 1;
            end
        end
    end
end

        
% keep calculating # mRNA in ROI.
for i = 1:size(rnas,1)
    if ~isempty(rnasNucList{i})
        nuc(rnasNucList{i},11) = nuc(rnasNucList{i},11)-1+1/length(rnasNucList{i});
    end
end
        
        
        
        
        
        
        
        

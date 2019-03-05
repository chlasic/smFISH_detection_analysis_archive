function [Nuc, thrNuc, aftNuc, tickn] = DetectNucleus(NucPla, iminfo, pix, tnran, sensi, tol, f, lic, thresForNuc)
% This function detects nuclei using circle detection method (hough
% transformation).

tickn = zeros(iminfo(4),1);
aftNuc = cell( iminfo(4), 1);
E = cell(iminfo(4),1);
F = cell(iminfo(4),1);
G = cell(iminfo(4),1);
H = cell(iminfo(4),1);
t1 = cell(iminfo(4),1);
temp = cell(iminfo(4),1);
f7 = cell(iminfo(4),1);
f8 = cell(iminfo(4),1);
mupr = 1 + thresForNuc;

for i = 1:iminfo(4)
    temp{i} = imadjust(NucPla{i,1}); 
    t1{i} = NucPla{i,1};

    tt = sort(t1{i}, 1);
    tt = sort(tt, 2);
    tickn(i) = 255/mean(mean(tt(end-2:end,end-2:end)))/ 3 ;

    if ismember(i,[1 2 3]) == 1
        tickn(i) = 1;
    end
    if ismember(i,[iminfo(4) iminfo(4)-1 iminfo(4)-2]) == 1
        tickn(i) = tickn(end-3);
    end
    if tickn(i) == 0
        tickn(i) = 1;
    end
    if tickn(i) > 4
        tickn(i) = tickn(i-1);
    end
end

parfor i = 1:iminfo(4)
    t2 = round(t1{i} * tickn(i)); 
    fint = adapthisteq(t1{i});
%     f1 = medfilt2(temp{i}, [round(pix*3.5) round(pix*3.5)]);  %%% old CHL setting.
    f1 = medfilt2(temp{i}, [round(pix*2) round(pix*2)]);  %%% new setting.

    f2 = imsharpen(f1, 'radius', pix*2, 'amount', pix);
       
    f3 = medfilt2(fint, [pix pix]);
    f3 = imsharpen(f3, 'radius', pix, 'amount', pix);
    f3 = imclose(f3, true(pix));        

    thrs = graythresh(f3) * mupr;
    if thrs > 1
        thrs = 1;
    end

    f5 = imbinarize(f3, thrs);
    f6 = imclose(f5, true(pix));
    f7n = bwareaopen(f6, pix*10);
    f7{i} = f7n;
    f8{i} = bwconvhull(f7n, 'objects');

    t3 = imcomplement(f7n);
    
    if isa(t2,'uint16')
        modr = uint16(t3);
    elseif isa(t2,'uint8')
        modr = uint8(t3);
    end
    t4 = mean(mean(t2 .* modr) .* tickn(i));

    aftNuc{i} = t2-t4;

    [A, B, C] = imfindcircles(f2, tnran, 'sensitivity', sensi);
    
    if isempty(A) == 0
        D = zeros(length(A(:,1)),1);
            for j=1:length(A(:,1))
                D(j,1) = GetIntensity(t2,A(j,:),B(j,:));
            end

        % remove dim (DAPI signal) nuclei
        cutInt = mean(D) * 0.6 ; %%%%% nuclear intensity below 'cutInt' will be removed.

        A(D < cutInt,:) = 0;
        B(D < cutInt,:) = 0;
        C(D < cutInt,:) = 0;
        D(D < cutInt,:) = 0;

        A = snip(A, '0');
        B = snip(B, '0');
        C = snip(C, '0');
        D = snip(D, '0');

        [E{i}, F{i}, G{i}, H{i}] = ...
            RemoveOverLap(A, B, C, D, tol, 4);

        E{i}(:,3) = 0;
    end
    
    fprintf('\nNucleus: %d(th)/total %d images,... %d(th)/ %d z-planes.', f, lic, i, iminfo(4));
end


NDNuc = [E F G H];
NDNuc(:,5) = {0};

parfor i = 1:iminfo(4)
    NDNuc{i,6} = i;
end

thrNuc = [f8 f7];


% variable NDNuc
% row: z-plane
% col 1: x,y coordinate of nuclear circles, col 2: radii of each circle
% col 3: circularity, col 4: total DAPI intensity of each circle, col 5: z-plane ID
NDNuc(any(cellfun(@isempty,NDNuc),2),:) = [];    

% fit nuclear circles into spheres.
% gather circles with the same (or close) centers (2nd method)
% 'cirs': col 1-2 | x, y coordinates; col 3 | z-plane; col 4 | diameter; col 5 | circularity; ; col 6 | total sig intensity
fgv = .55  / mean(iminfo(5:6,1)); %%%%% 0.6 um distance between two centers of circle is accepted to fit one sphere.
id = 0; % sphere ID number
bfr = round(  0.7  /iminfo(7)); %%%%% will use circles in 1 um z-range to fit sphere.
cirs = cell(9999,1);



for i = 1: length(NDNuc(:,1))-bfr % each plane
    for j = 1:length(NDNuc{i,1}(:,1)) % each nucleus
        if NDNuc{i,1}(j,3) == 0
            id = id + 1;
            NDNuc{i,1}(j,3) = id;
%             counter(id,2) = 1;
            idNow = id;

            % saves circles in a sphere in one cells.
            cirs{idNow,1} = [NDNuc{i,1}(j,1:2) NDNuc{i,6} NDNuc{i,2}(j,1) NDNuc{i,3}(j,1) NDNuc{i,4}(j,1)];

        else
            idNow = NDNuc{i,1}(j,3);
        end

        x = NDNuc{i,1}(j,1); y = NDNuc{i,1}(j,2);

        tmp = cirs{idNow,1}(:,3);
        loc = i+1:i+bfr;
        chan = cell2mat(NDNuc(loc,6));
        loc(ismember(chan,tmp)) = []; 
        for k = loc % plane for scan
            for l = 1:length(NDNuc{k,1}(:,1)) % nucleus for scan
                if length(cirs{idNow,1}(:,1)) > 4 &&...
                        NDNuc{k,2}(l,1) > cirs{idNow,1}(end,4) && NDNuc{k,4}(l,1) > cirs{idNow,1}(end,6)*1.1
                elseif NDNuc{k,2}(l,1) / cirs{idNow,1}(end,4) < 1.25 && NDNuc{k,2}(l,1) / cirs{idNow,1}(end,4) > 0.8
                    a = NDNuc{k,1}(l,1); b = NDNuc{k,1}(l,2);
                    if cirs{idNow,1}(end,5) > 0.1 && NDNuc{k,3}(l,1) > 0.1
                        if sqrt((x-a)^2 + (y-b)^2) < fgv 
                            NDNuc{k,1}(l,3) = NDNuc{i,1}(j,3);
                            ntmp = [NDNuc{k,1}(l,1:2) NDNuc{k,6} NDNuc{k,2}(l,1) NDNuc{k,3}(l,1) NDNuc{k,4}(l,1)];
                            cirs{idNow,1} = [cirs{idNow,1}; ntmp ];
                        end
                    else
                        if sqrt((x-a)^2 + (y-b)^2) < fgv*3 
                            NDNuc{k,1}(l,3) = NDNuc{i,1}(j,3);
                            ntmp = [NDNuc{k,1}(l,1:2) NDNuc{k,6} NDNuc{k,2}(l,1) NDNuc{k,3}(l,1) NDNuc{k,4}(l,1)];
                            cirs{idNow,1} = [cirs{idNow,1}; ntmp ];
                        end
                    end
                end
            end
        end

    end
end

if id == 0
    id = 1;
end
cirs(id:end,:) = [];




% calculate size of nuclear circles (detected) and find its center.
% 'cirs_fin': | x-coor (in pixel) | y | nth z-plane | radius of circle |
%               | circularity (0-1) | total DAPI signal in circle |
cir_size = cellfun(@numel, cirs)/6;
idx = cir_size > 2; %%%%% only circle seen more than 3 slides are considered.
cirs_fin = cirs;
cirs_fin = cirs_fin(idx,1);

% generate a sphere from detected circles on the same nucleus.
% 'Nuc' : | x-coor | y | z | radius | ave. circularity |
%           | std. circularity | total DAPI sig. | 
Nuc = zeros (length(cirs_fin(:,1)),4);

for i= 1:length(cirs_fin(:,1))
    len = length(cirs_fin{i,1}(:,1));
    circoor = zeros(len*4,3);
    for j = 1:len
        circoor(4*j-3:4*j,:) = [cirs_fin{i,1}(j,1)+cirs_fin{i,1}(j,4) cirs_fin{i,1}(j,2) cirs_fin{i,1}(j,3)*iminfo(7)/iminfo(5);...
            cirs_fin{i,1}(j,1) cirs_fin{i,1}(j,2)+cirs_fin{i,1}(j,4) cirs_fin{i,1}(j,3)*iminfo(7)/iminfo(5);...
            cirs_fin{i,1}(j,1)-cirs_fin{i,1}(j,4) cirs_fin{i,1}(j,2) cirs_fin{i,1}(j,3)*iminfo(7)/iminfo(5);...
            cirs_fin{i,1}(j,1) cirs_fin{i,1}(j,2)-cirs_fin{i,1}(j,4) cirs_fin{i,1}(j,3)*iminfo(7)/iminfo(5)];
    end

    % generate coordinates of dots on the circle to use 'sphereFit.'
    [Nuc(i,1:3),Nuc(i,4)] = sphereFit(circoor);
    Nuc(i,5:7) = [ mean(cirs_fin{i,1}(:,5)) std(cirs_fin{i,1}(:,5)) sum(cirs_fin{i,1}(:,6)) ];

end
    
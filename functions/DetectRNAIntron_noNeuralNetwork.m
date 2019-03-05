function [blobrna, blrna, mom] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, lic, g, thresForIntron)
% This function detects intronic RNA spots (transcription sites). The
% output is RNA coordinates, counts and intensity


sizel = iminfo(4);
fint = cell(sizel,1);  
t1 = cell(sizel,1); % cropped RNA images only in the nuclei.
mofm = zeros(sizel,1);
thrrna = cell( sizel, 1);
sen = strel('disk', pixr*15); 
se = strel('disk', pixr*2); 
mupr = thresForIntron+0.5; % multiplier to define the background level (becomes higher when bg is high) 
                        
% make general threshold
parfor i = 1:sizel
    t1{i} = rPla{i,1};
    
    % black RNA images not overlapping with nuclei
    if ismember(i, [1 2 3]) == 1
        tb = cat(3,thrNuc{1:i+3,2})*1.5;
    elseif ismember(i, [iminfo(4)-2 iminfo(4)-1 iminfo(4)]) == 1
        tb = cat(3,thrNuc{i-3:end,2})*1.5;
    else
        tb = cat(3,thrNuc{i-3:i+3,2});
    end
    tb = max(tb,[],3);
    tx = imdilate(tb,sen);
    
    if isa(t1{i},'uint16')
        modr = uint16(tx);
    elseif isa(t1{i},'uint8')
        modr = uint8(tx);
    end
    t1{i} = t1{i} .* modr;
    
    % make a binary image where background signal is.
    thrs = graythresh(rPla{i,1})/3;
    tmask = im2bw(t1{i}, thrs * 1 );
    
    % calculate background level only inside the gonad using gonadal mask.
    % calculate background for individual z planes and save in 'mofm'
    tAll = t1{i}(tmask > 0);
    mofm(i) = mean(tAll)*1.5;  
end

% remove NaN from 'mofm' to calculate 'mmofm'
mofm(isnan(mofm)) = 1;
mmofm = mean(mofm);
mom = mmofm;

mofm = mofm * mupr;
mmofm = mmofm * mupr;



%%% threshold images and detect RNA spots
parfor i = 1:sizel
    temp = t1{i};
    fint{i} = (rPla{i,1} - mofm(i)/mupr) / mofm(i)*mupr * mmofm;

    temp(temp < mofm(i) ) = 0;

    t2 = temp; 
    
    if mom >= 2.5 
        loopn = 2;
    else
        loopn = 1;
    end
    
%     loopn = 1;

    for j = 1:loopn
        t2 = medfilt2(t2, [pixr pixr]);
        tt = reshape(t2,[],1);
        tt = sort(tt);
        ty = tt(tt ~= 0);
        tmofm = mean(mean(ty(20:end)));
        t2(t2 < tmofm) = 0;
    end

    % correct pixel registration
%     tpix = t2(end-1:end,:);
%     t2 = t2(1:end-2,:);
%     t2 = [tpix; t2];
%     tpix = t2(:,end-1:end);
%     t2 = t2(:,1:end-2);
%     t2 = [tpix t2];
    
    tpix = t2(end,:);
    t2 = t2(1:end-1,:);
    t2 = [tpix; t2];
    tpix = t2(:,end);
    t2 = t2(:,1:end-1);
    t2 = [tpix t2];

% ------------- peak detection method (1/3) ------------------
    [A, B] = FastPeakFind(t1{i}, mofm(i)*mupr*0.8);
    A = [A(1:2:end) A(2:2:end)];
        

    % remove peaks too close
    A(:,3) = 0;
    for j = 1:length(A(:,1))
        for k = j+1:length(A(:,1))
            dist = sqrt((A(k,1)-A(j,1))^2 + (A(k,2)-A(j,2))^2);
            if dist < pixr*2
                A(j,3) = 999;
                B(A(j,2),A(j,1)) = 0;
            end
        end
    end
    A(A(:,3) == 999,:) = [];
    A(:,3) = [];
    
    A(:,3) = A(:,1);
    A(:,1) = [];
    
%=========== visual detected RNAs =======
% figure,imshow(10*t1{i});
% hold on
% plot(A(:,2), A(:,1), 'r+');
% pause;
% close all
% end
% ---------------------------------------


% ===========%%% bwconncomp method (2/3)==============
    thrs = graythresh(fint{i});
    if thrs > 1
        thrs = 1;
    elseif thrs == 0
        thrs = 0.01;
    end

    tcou = 0;
    gogo = 1;
    
    if mom >= 2.5 
        m2 = 2;
    else
        m2 = 1;
    end
    
    while gogo == 1 && tcou < 100
        t3 = im2bw(t2, thrs * .5 );
        t5 = bwareaopen(t3, ceil(pixr*m2));
        t6 = t5 & ~bwareaopen(t5, ceil(pixr*30));
        tgo = bwconncomp(t6, 8);
        if tgo.NumObjects > 40
            thrs = thrs*1.1;
            tcou = tcou + 1;
        else
            gogo = 0;
        end
    end
    t6 = imdilate(t6, se);
    %---------------------------------------------------------  
    
    %%% cross validate detected spots (bwconncomp vs. peakfinder)
    % Remove bwconncomp-detected dots that are not found by peakfinder.
    if ~isempty(A) && sum(sum(A)) ~= 0
        A(:,3) = 0;
        for j=1:size(A,1)
            if t6(A(j,1),A(j,2)) ~= 0
                A(j,3) = 999; % mark valid peaks
            else
                B(A(j,1),A(j,2)) = 0;
            end
        end
        A(A(:,3) == 0,:) = [];
        A(:,3) = [];
    end


    % make a binary image where background signal is.
%     thrs = graythresh(rPla{i,1})/3;
%     tmask = im2bw(t1{i}, thrs * 1 );
%     
%     se2 = strel('disk', pixr*4); 
%     se3 = strel('disk', pixr*15); 
% 
%     tm2 = imdilate(tmask,se2);
%        
%     tm3 = imerode(tm2,se3);
%     
% 
%     % get intensity of spots and background to get the s:n ratio.
%     sizC = pixr*2-1; % radius of the center region to get intensity
%     cal = zeros(length(A(:,1)), 4);
%     if ~isempty(A) && sum(sum(A)) ~= 0
%         for wsz = [2 3 4]
%             for j=1:length(A(:,1))
%                 % 'A' col: | x-coor | y-coor | center_intensity | backgd2 int | bkgd3 int | bkgd4 int |.
%                 sizB = pixr*2*wsz; 
%                 xb = 1:sizB*2+1;
%                 
%                 [x, y] = meshgrid(xb,xb);
%                 Rb = sqrt((x-sizB-1).^2 + (y-sizB-1).^2); % center circle
%                 
%                 Rc = Rb;
%                 Rc(Rc>sizC) = 0;
%                 Rc(sizB+1,sizB+1) = 1;
%                 Rb(Rb>sizB) = 0;
%                 Rb(Rb<sizC+2) = 0;
%                 
%                 LB = A(j,1) - sizB; % left boundary of the ROI for scanning
%                 if LB < 1
%                     Rb(:,1:1-LB)=[];
%                     Rc(:,1:1-LB)=[];
%                     LB = 1;
%                 end
%                 
%                 RB = A(j,1) + sizB;
%                 if RB > iminfo(3)
%                     Rb(:,end-RB+iminfo(3)+1:end)=[];
%                     Rc(:,end-RB+iminfo(3)+1:end)=[];
%                     RB = iminfo(3);
%                 end
%                     
%                 UB = A(j,2) - sizB;
%                 if UB < 1
%                     Rb(1:1-UB,:)=[];
%                     Rc(1:1-UB,:)=[];
%                     UB = 1;
%                 end
%                 
%                 DB = A(j,2) + sizB;
%                 if DB > iminfo(2)
%                     Rb(end-DB+iminfo(2)+1:end,:)=[];
%                     Rc(end-DB+iminfo(2)+1:end,:)=[];
%                     DB = iminfo(2);
%                 end
%                 
%                 cutcen = rPla{i,1}(LB:RB,UB:DB);
%                 centr = cutcen(Rc > 0);
%                 
%                 cutbck = rPla{i,1}(LB:RB,UB:DB);
%                 surrd = cutbck(Rb > 0);
%     
% 
% 
%                 % subtract detected spots that are dim (not significant spots).
%                 [~,cal(j,1)] = kstest2(centr, surrd);
%                 cal(j,2) = mean(centr)/mean(surrd); %ratio to surrounding bg.
%                 cal(j,3) = mean(centr)/mofm(i); % ratio to mean background (fixed)
%                 cal(j,4) =  j;
%                 caller = 1;
% 
% %                 if tm3(A(j,1),A(j,2)) == 0
% %                     B(A(j,1),A(j,2)) = 0;
% %                     caller = 0;
% %                 else
%                 
%                 if wsz == 2
%                     if cal(j,1) > 1*10^-(mupr * 0) || cal(j,2) < mupr * .8 %|| cal(j,3) < mupr * 0.5 %|| A(j,2) > iminfo(2)/2.5
%                         B(A(j,1),A(j,2)) = 0;
%                         caller = 0;
%                     end
%                 elseif wsz == 3
%                     if cal(j,1) > 1*10^-(mupr * 0) || cal(j,2) < mupr * .8
%                         B(A(j,1),A(j,2)) = 0;
%                         caller = 0;
%                     end
%                 elseif wsz == 4
%                     if cal(j,1) > 5*10^-(mupr * 0) || cal(j,2) < mupr * .8
%                         B(A(j,1),A(j,2)) = 0;
%                         caller = 0;
%                     end
%                 end
%                 
% %                 end
% 
%                 % =========== show results ========================
% %                 fprintf('\n%dth-> p: [[%d]] %d, ratio: %d, ratio(pos back): %d', j, caller, cal(j,1), cal(j,2), cal(j,3));
% %                 
% %                 if j==length(A(:,1))
% %                     fprintf('\n   ^^^ < %dXradius applied >',wsz);
% %                     fprintf('\n');
% %                 end
%                     %---------- save results for review -----------
% 
% 
%                 % ==================================================
%             end
%         end   
%     end
    
    tB = imdilate(B, strel('disk',pixr*2+1));
    
    thrrna{i} = tB;
    fprintf('\nRNA%d: %d(th)/total %d images,... %d(th)/ %d z-planes.', g, f, lic, i, sizel);
%     figure,imshow(t7) % use it with for loop
% imshow(tB);
% hold on
% plot(A(:,1), A(:,2), 'r+');
% pause
end



% check result (binary images with only meaningful spots)
% figure
% for i = 1:sizel
%     imshow(thrrna{i})
%     fprintf('\n%d -th stack',i);
%     pause
% end




%%% 3D reconstitution using bwconncomp.
% blrna: |Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
%       | 6th: matching nuc #(row # in 'nuc')
tfrna = cat(3,thrrna{:});
temp = bwconncomp(tfrna, 26);
conrna = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');
blrna = struct2cell(conrna)';

    
% calculates intensity of detected blobs from original images.
iint = zeros(1,length(blrna(:,1)));
for i=1:length(blrna(:,1))
    iint(i) = 0;
    for j=1:blrna{i,3}(6)
        cutimg = fint{round(blrna{i,3}(3)),1};
        cutimg = cutimg(round(blrna{i,3}(2)):round(blrna{i,3}(2))+blrna{i,3}(5)-1, round(blrna{i,3}(1)): round(blrna{i,3}(1))+blrna{i,3}(4)-1);

        iint(i) = iint(i) + sum(sum(blrna{i,4}(:,:,j) .* double(cutimg)));
    end
end

blrna = [blrna num2cell(iint')];
% fprintf('\nsize: %d %d\n', size(blrna,1), size(blrna,2));    



if  isempty(blrna) == 0
    if isempty(blrna{1,1}) == 0
        % 'blobrna': | rna ID | x-size | y-size | z-size | 5th: total # pixel | 
        %           | total intensity | x centroid | y centroid | z centroid | # PTS in the blob.
        blobrna = zeros(length(blrna(:,1)),1);
        blobrna(:,1) = 1:length(blrna(:,1));
        temp = cell2mat(blrna(:,3));
        blobrna(:,2:4) = temp(:,4:6);
        temp = cell2mat(blrna(:,2));
        blobrna(:,7:9) = temp(:,1:3);

        blobrna(:,6) = cell2mat(blrna(:,5));
        blobrna(:,5) = cell2mat(blrna(:,1));

            
        %%% remove false-positively detected blobs
        if isempty(blobrna) == 0
            % get mean of obj. size (pixels) and intensity (normalized: ind. int. - mean)
            mint = mean(blobrna(:,6));
            mpix = mean(blobrna(:,5));

            % take out objects detected less than 3 z-planes
            % take out objects dim/small.
            temp = blobrna(:,4) < 2 & blobrna(:,6) < mint*0.1 ;
            blobrna(temp,:) = [];
            
            
            if isempty(blobrna) == 0
                temp = blobrna(:,4) < 2 & blobrna(:,5) < mpix*0.1 ;
                blobrna(temp,:) = [];
                
                if isempty(blobrna) == 0
                    % take out < 3x3 obj.
                    temp = blobrna(:,2) < pixr+1 & blobrna(:,3) < pixr+1;
                    blobrna(temp ,:) = [];

%                     if isempty(blobrna) == 0
%                         % take out obj with total intensity < 100 a.u. (possibly just background or bleed-through)
%                         temp = blobrna(:,6) < mmofm*10;
%                         blobrna(temp,:) = [];
%                     end

%                     if isempty(blobrna) == 0 && f == 12
%                         % take out obj on the edge
%                         temp = (blobrna(:,7) < 486 & blobrna(:,8) < 174);
% %                         temp = (blobrna(:,7) < 130 & blobrna(:,8) < 130) | ...
% %                         (blobrna(:,7) < 255 & blobrna(:,8) > iminfo(3)-130) | ...
% %                         (blobrna(:,7) > iminfo(2) - 130 & blobrna(:,8) < 130) | ...
% %                         (blobrna(:,7) > iminfo(2) - 130 & blobrna(:,8) > iminfo(3)-130);
% 
%                         blobrna(temp,:) = [];
%                     end
                end
            end
        end
    else
        blobrna = [];
    end
else
    blobrna = [];
end












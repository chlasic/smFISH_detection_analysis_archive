function [cyl_axis, outlc, t2] = DefOutlineNTS(mip, pix, iminfo)
% This function finds cell outlines based on input phase images and outputs
% coordinates and radii (to the outline) of cell center line. Input
% variables are following. PhaPla: phase image stack (taken out from
% confocal images), pix: smallest pixel size control (automatically 
% determined), iminfo: image info imported from metafile (automatically
% done).


% apply image filters for crisper images
tm1 = imadjust(mip);
t1 = medfilt2(tm1, [pix*2 pix*2]);
t1(t1<mean(mean(t1))/3) = 0;
t1(t1>0) = 1;
% t2 = imcomplement(t1); % invert image
t2 = logical(t1);

outcoor = bwboundaries(t2);
outsiz = cell2mat(cellfun(@length, outcoor, 'uni', 0));
outsiz = outsiz(:,1);
% outlc contains x,y coor of germline outline.
outlc = cell2mat(outcoor(outsiz == max(outsiz)));
outlc(:,[1,2]) = outlc(:,[2,1]);

% reconsititute as 3D cylinder
outlc(:,3) = round(iminfo(4)/2);

%%% gets cylinderical axis coordinates and r to draw the cell on the plot.
%%% "cline" is a sorted Cell outline coordinates.
%%% (1) x-coor | (2) y-coor | (3) z-coor | (4) Nb mask points on that x |
%%%
%%% "cyl_axis" contains: (unit: pixels)
%%% (1) x-coor | (2) y-coor | (3) z-coor | (4) radius of the cell on that point.| 
cline = sortrows(outlc, 1);
cline(:,4) = 1;
cline = sortrows(cline, [1 2]);

for i=2:length(cline(:,1))
    if cline(i,1) == cline(i-1,1)
        cline(i,4) = cline(i,4) + cline(i-1,4);
    end
end
cline = sortrows(cline, [1 -4]);

lent = cline(end,1)-cline(1,1)+1; 
jr = 1;
cyl_axis= zeros(lent,4);
for i=1: lent
    cyl_axis(i,1:4) = [cline(jr,1) (cline(jr,2)+cline(jr+cline(jr,4)-1,2))/2 cline(jr,3) (cline(jr,2)-cline(jr+cline(jr,4)-1,2)/2)/2];
    jr = jr + cline(jr,4);
end
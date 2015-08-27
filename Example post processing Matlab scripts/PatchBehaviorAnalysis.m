scrsz = get(0,'screensize');
clr1 = [0.2 0.1 0.8];
clr2 = [0.8 0.1 0.1];
clr3 = [0.1 0.2 0.3];
clr4 = [0.2 0.8 0.5];
clr5 = [0.6 0.1 0.6];
clr6 = [0.9 165/255 0];

frameRate = 5;

%%
nosetip = [];
PatchResidency = zeros(3,size(RT_L,2))*NaN;
NoseToPatchDistances = zeros(NumPatches,size(RT_L,2))*NaN;
t_WL =nanmean(WormLengthRT);
Closeness = 0.05; % What distance from patch edge to be counted as patch residency in fraction of worm length

for i = 1:size(RT_L,2)
    nosetip(i,:) = nosetailRT{i}(1,:);
    if isfinite(nosetip(i,1))
        t1 = nosetip(i,:);
        t_flag = 0;
        t_dist_Coll = [];
        for j = 1:NumPatches
            t_dist = pdist2(t1,PatchPixels{j});
            t_dist = min(t_dist);
            
            NoseToPatchDistances(j,i) = t_dist;
            
            t_dist_Coll(j,:) = [j t_dist];
            if t_dist < t_WL*Closeness
                
                PatchResidency(1,i) = 1;
                PatchResidency(2,i) = j;
                PatchResidency(3,i) = t_dist;
                t_flag = 1;
            end
        end
        if t_flag == 0
            PatchResidency(1,i) = 0;
            PatchResidency(2,i) = find(t_dist_Coll(:,2) == min(t_dist_Coll(:,2)));
            PatchResidency(3,i) = min(t_dist_Coll(:,2));
        end
    end
end

%%
CoM_Vel = zeros(1,size(RT_L,2))*NaN;
t_step = 6; %frames
for i = floor((t_step)/2) + 1:size(RT_L,2)-floor((t_step)/2)
    CoM_Vel(i) = (pdist2(CoM(i-floor((t_step)/2),:),CoM(i+(floor((t_step)/2)),:)))/(t_step/frameRate);
end
%%
PatchBouts = {};
indx = isfinite(PatchResidency(1,:));
PatchResidencyTemp = PatchResidency(:,indx);
    
timeTemp = [1:size(PatchResidency,2)];
timeTemp = timeTemp(indx);
    
InPatchEdges = find(diff(PatchResidencyTemp(1,:)));
InPatchTemp = [];
if ~isempty(InPatchEdges)
    if PatchResidencyTemp(1,InPatchEdges(1)) == 1
        InPatchTemp = [InPatchTemp; [1 timeTemp(InPatchEdges(1))]];
    end
end
for k = 1:length(InPatchEdges)-1
    if PatchResidencyTemp(1,InPatchEdges(k)+1) == 1
        InPatchTemp = [InPatchTemp; [timeTemp(InPatchEdges(k)+1) timeTemp(InPatchEdges(k+1))]];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Filter out short 2 second out of patch excursions
InPatchTempF = []; k = 0;
if ~isempty(InPatchTemp)
    k = k + 1;
    InPatchTempF(k,:) = InPatchTemp(1,:);
    for i = 2:size(InPatchTemp,1)
        x = InPatchTemp(i,1)-InPatchTemp(i-1,2);
        if x > 0*frameRate
            k = k + 1;
            InPatchTempF(k,:) = InPatchTemp(i,:);
        else
            InPatchTempF(k,2) = InPatchTemp(i,2);
        end
    end
end
%%%Filter out short 2 second inPatch 
InPatchF = []; k = 0;
for i = 1:size(InPatchTempF,1)
    if diff(InPatchTempF(i,:)) > 0*frameRate
        k = k + 1;
        InPatchF(k,:) = InPatchTempF(i,:);
    end
end
InPatchCross = InPatchF;
%%
%%%%%%%%%%%%%%
CoM_Vel_Roam = CoM_Vel;
for i = 1:size(InPatchF,1)
        CoM_Vel_Roam(InPatchF(i,1):InPatchF(i,2)) = NaN;
end
%%
% 
InPatchF2 = []; k = 0;
if ~isempty(InPatchF)
    k = k + 1;
    InPatchF2(k,:) = InPatchF(1,:);
    for i = 2:size(InPatchF,1)
        x = InPatchF(i,1)-InPatchF(i-1,2);
        
        range1 = InPatchF(i,1):InPatchF(i,2);
        range2 = InPatchF(i-1,2):InPatchF(i,1);
        
        t1 = nanmean(PatchResidency(2,range1)); 
        t2 = nanmean(PatchResidency(2,range2));
        
        t3 = max(PatchResidency(3,range2));
        if t3 > t_WL * 0.1 || t1 ~= t2
            k = k + 1;
            InPatchF2(k,:) = InPatchF(i,:);
        else
            InPatchF2(k,2) = InPatchF(i,2);
        end
    end
end
InPatchF = InPatchF2;

InPatchF2 = []; k = 0;
for i = 1:size(InPatchF,1)
    Dist_Patch_t = PatchResidency(3,InPatchF(i,1):InPatchF(i,2));
    Dist_Patch_t = Dist_Patch_t(isfinite(Dist_Patch_t));
    FracInPatch = sum(Dist_Patch_t<2)/length(Dist_Patch_t);
    t1 = nanmean(PatchResidency(2,InPatchF(i,1):InPatchF(i,2)));
    %dwl_vel = (nanmean(CoM_Vel_Roam)/3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALPHA = 0.15;
    tmp = sort(CoM_Vel_Roam(isfinite(CoM_Vel_Roam)));
    tmpL = length(tmp);
    dwl_vel = tmp(round(ALPHA*tmpL));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(t1,1) ~= 0
        continue
    end
    DurFilt = ((PatchMajAxis(t1)+ 2*(t_WL*Closeness)) / dwl_vel) * frameRate; % filter out Short periods
    
    if FracInPatch > 0.5
        if diff(InPatchF(i,:)) < DurFilt
        else
            k = k + 1;
            InPatchF2(k,:) = InPatchF(i,:);
        end
    end
end
InPatchF = InPatchF2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%FIlter out Out of Patch Excursions that Don't make it further then 20% of WL from Patch %%%%%%%%%%%%%%
InPatchF2 = []; k = 0;
if ~isempty(InPatchF)
    k = k + 1;
    InPatchF2(k,:) = InPatchF(1,:);
    for i = 2:size(InPatchF,1)
        x = InPatchF(i,1)-InPatchF(i-1,2);
        
        range1 = InPatchF(i,1):InPatchF(i,2);
        range2 = InPatchF(i-1,2):InPatchF(i,1);
        
        t1 = nanmean(PatchResidency(2,range1)); 
        t2 = nanmean(PatchResidency(2,range2));
        
        t3 = max(PatchResidency(3,range2));
        if t3 > t_WL * 0.2 || t1 ~= t2
            k = k + 1;
            InPatchF2(k,:) = InPatchF(i,:);
        else
            InPatchF2(k,2) = InPatchF(i,2);
        end
    end
end
%%%%%%%%%%%%%%%
InPatch_cross = []; k = 0;
for i = 1:size(InPatchCross,1)
    if diff(InPatchCross(i,:)) > frameRate
        range = InPatchCross(i,1):InPatchCross(i,2);
        t1 = nanmean(PatchResidency(2,range));
        if mod(t1,1) ~= 0
            continue;
        end
        k = k + 1;
        InPatch_cross(k,:) = InPatchCross(i,:);
    end
end

%%%%%%%%%%%%%%%

InPatchFilt = InPatchF2;

%%%%%%%%%%%%%%
InPatch_L = zeros(1,size(RT_L,2));
OutPatch_L = ones(1,size(RT_L,2));

for i = 1:size(InPatchFilt,1)
    InPatch_L(InPatchFilt(i,1):InPatchFilt(i,2)) = 1;
    OutPatch_L(InPatchFilt(i,1):InPatchFilt(i,2)) = 0;
end
InPatch_L = logical(InPatch_L);
OutPatch_L = logical(OutPatch_L);

OutPatchEdges = find(diff(OutPatch_L));

OutPatchFilt = []; k = 0;
if OutPatch_L(OutPatchEdges(1)) == 1
    OutPatchFilt = [OutPatchFilt; [1 OutPatchEdges(1)]];
end
for i = 1:size(OutPatchEdges,2)-1
    if OutPatch_L(OutPatchEdges(i)+1) == 1
        OutPatchFilt = [OutPatchFilt; [OutPatchEdges(i)+1 OutPatchEdges(i+1)]];
    end
end

%%%%%%%%%%%%%%
PatchBouts.InPatchIndx = InPatchFilt;
PatchBouts.OutPatchIndx = OutPatchFilt;
if InPatchFilt(1,1) < OutPatchFilt(1,1)
    PatchBouts.First = 'In';
else
    PatchBouts.First = 'Out';
end
PatchBouts.InPatch_L = InPatch_L;
PatchBouts.OutPatch_L = OutPatch_L;
%%%%%%%%%%%%%%%%%%
PatchBouts.PatchCross = InPatch_cross;
%%
InPatch_ID = [];
InPatch_Dist = [];
InPatch_PercInFood = [];
t = [];t2 = zeros(1,size(RT_L,2))* NaN;

InPatch_Time = [];%zeros(1,size(RT_L,2))*NaN;
for i = 1:size(PatchBouts.InPatchIndx,1)
    range = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
    t2(range) = (PatchResidency(3,range));
    
    InPatch_ID(i) = nanmean(PatchResidency(2,range));
    InPatch_Dist(i) = nanmean(PatchResidency(3,range));
    
    t = PatchResidency(3,range); t = t(isfinite(t));
    t2 = sum(t<2);
    
    InPatch_PercInFood(i) = t2/length(t) *100;
    
    t(i,1) = mean(range);
    t(i,2) = max(PatchResidency(3,range));
    t(i,3) = nanmean(PatchResidency(2,range));
    
    InPatch_Time(i,:) = [(mean(range)) length(range)];
    
end

PatchBouts.InPatchID = InPatch_ID;
PatchBouts.InPatchDur = diff(PatchBouts.InPatchIndx')/frameRate;
PatchBouts.InPatchDist = InPatch_Dist;
PatchBouts.InPatchPercInFood = InPatch_PercInFood;
PatchBouts.OutPatchDur = diff(PatchBouts.OutPatchIndx')/frameRate;
%%%%%%%%%%%%%%%

tt = [];
PatchBouts.NewPatchIndx(1) = 1;
for i = 2:size(PatchBouts.InPatchIndx,1)
    range = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
%     t(i) = round(nanmean(PatchResidency(2,range)));
    tt = PatchResidency(2,range);
    tt = tt(isfinite(tt));
    tt_i = (tt==PatchBouts.InPatchID(i-1));
    tt_p = sum(tt_i)/length(tt_i);
    if tt_p > 0.5 && PatchBouts.InPatchID(i-1) == PatchBouts.InPatchID(i)
        PatchBouts.NewPatchIndx(i) = 0;
    else
        PatchBouts.NewPatchIndx(i) = 1;
    end
end


%%%%%%%%%%%%%%%%
CrossPatch_ID = [];
CrossPatch_Vel = [];
for i = 1:size(PatchBouts.PatchCross,1)
    range = PatchBouts.PatchCross(i,1):PatchBouts.PatchCross(i,2);
    CrossPatch_ID(i) = nanmean(PatchResidency(2,range));
    CrossPatch_Vel(i) = nanmean(CoM_Vel(range));
end
PatchBouts.CrossPatchID = CrossPatch_ID;
PatchBouts.CrossPatchDur = diff(PatchBouts.PatchCross')/frameRate;
PatchBouts.CrossPatchVel = CrossPatch_Vel;

C = colormap(jet(NumPatches)); close;
figure;hold on;
plot(InPatch_Time(:,1),InPatch_Time(:,2)/frameRate,'k.')
for i = 1:size(t,1)
    text(InPatch_Time(i,1)+5,InPatch_Time(i,2)/frameRate+5,num2str(PatchBouts.InPatchID(i)),'color',C(PatchBouts.InPatchID(i),:))
end
set(gcf,'position',scrsz);
saveas(gcf,[Worm  ' InPatch_Duration.pdf']);
close
%%
InPatch_Residency = cell(NumPatches,1);
for i = 1:NumPatches
    InPatch_Residency{i}.EncounterNum = 0;
    InPatch_Residency{i}.Num = 0;
    InPatch_Residency{i}.CumDurs = [0];
end

for i = 1:size(PatchBouts.InPatchIndx,1)
    InPatch_Residency{PatchBouts.InPatchID(i)}.Num = InPatch_Residency{PatchBouts.InPatchID(i)}.Num + 1;
    InPatch_Residency{PatchBouts.InPatchID(i)}.Frames(InPatch_Residency{PatchBouts.InPatchID(i)}.Num,:) = PatchBouts.InPatchIndx(i,:);
    InPatch_Residency{PatchBouts.InPatchID(i)}.Durs(InPatch_Residency{PatchBouts.InPatchID(i)}.Num) = PatchBouts.InPatchDur(i);
    InPatch_Residency{PatchBouts.InPatchID(i)}.CumDurs(InPatch_Residency{PatchBouts.InPatchID(i)}.Num) = (InPatch_Residency{PatchBouts.InPatchID(i)}.CumDurs(end))+PatchBouts.InPatchDur(i);

    InPatch_Residency{PatchBouts.InPatchID(i)}.NewPatch(InPatch_Residency{PatchBouts.InPatchID(i)}.Num) = PatchBouts.NewPatchIndx(i);
end


t_max = [];
for i = 1:NumPatches
    t_max(i) = InPatch_Residency{i}.Num;
end
t_max = max(t_max);
% C = colormap(jet(t_max)); close;

temp_Res = zeros(NumPatches,t_max);
for i = 1:NumPatches
    if InPatch_Residency{i}.Num ~=0
    for j = 1:size(InPatch_Residency{i}.Durs,2)
        temp_Res(i,j) = InPatch_Residency{i}.Durs(j);
    end
    end
end

figure;
bar(temp_Res,'stacked')
set(gcf,'position',scrsz);
xlabel('Patch ID')
ylabel('Duration (sec)')
saveas(gcf,[Worm  ' InPatch_Residency.pdf']);
close

% figure;bar(fliplr(sort(sum(temp_Res'))))

%%
for i = 1:size(PatchBouts.PatchCross,1)
    InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterNum = InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterNum + 1;
    InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterFrames(InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterNum,:) = PatchBouts.PatchCross(i,:);
    InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterDurs(InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterNum) = PatchBouts.CrossPatchDur(i);
    InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterVels(InPatch_Residency{PatchBouts.CrossPatchID(i)}.EncounterNum) = PatchBouts.CrossPatchVel(i);
end

for k = 1:size(InPatch_Residency,1)
    tmp_patch_pixels = PatchPixels{k};
    FeedingArea = [];
    for i = 1:InPatch_Residency{k}.Num
        range = InPatch_Residency{k}.Frames(i,1):InPatch_Residency{k}.Frames(i,2);
        p_bout_pix = round([nosetip(range,1) nosetip(range,2)]);
        
        p_bout_img = (zeros(2040,2050));
        
        
        p_pix_indx = p_bout_pix(:,2)+((p_bout_pix(:,1)-1)*size(p_bout_img,1)); % linear index into bw
        p_pix_indx = p_pix_indx(isfinite(p_pix_indx));
        p_bout_img(p_pix_indx) = 1;
        
        ttt = bwmorph(p_bout_img,'bridge');
        ttt = bwmorph(ttt,'diag');
        ttt = bwmorph(ttt,'dilate');
        ttt = bwmorph(ttt,'dilate');
        ttt = bwmorph(ttt,'dilate');
        ttt = bwmorph(ttt,'erode');
        ttt = bwmorph(ttt,'erode');
        ttt = bwmorph(ttt,'majority');
        ttt = bwmorph(ttt,'close');
        
        [L,NUM] = bwlabel(ttt);
        STATS = regionprops(L,{'Area','PixelList'});
        indx_2 = find([STATS.Area] == max([STATS.Area])); indx_2 = indx_2(1);
        
        InPatch_Residency{k}.FeedingAreas{i} = STATS(indx_2).PixelList;
        match_pix = ismember(tmp_patch_pixels,STATS(indx_2).PixelList,'rows');
        InPatch_Residency{k}.PercEaten(i) = (sum(match_pix)/length(match_pix))*100;
       
        FeedingArea = [FeedingArea; STATS(indx_2).PixelList];
        
        match_pix = ismember(tmp_patch_pixels,FeedingArea,'rows');
        InPatch_Residency{k}.CumPercEaten(i) = (sum(match_pix)/length(match_pix))*100;
    end
    InPatch_Residency{k}.CumFeedingArea = FeedingArea;
    
    
    if ~isempty(FeedingArea)
        InPatch_Residency{k}.TotFeedingArea = FeedingArea;
        match_pix = ismember(tmp_patch_pixels,FeedingArea,'rows');
        InPatch_Residency{k}.TotPercEaten = (sum(match_pix)/length(match_pix))*100;
    else
        InPatch_Residency{k}.TotFeedingArea = 0;
        InPatch_Residency{k}.TotPercEaten = 0;
    end
    if InPatch_Residency{k}.EncounterNum ~= 0
        enc_bout_pix = [];
        for i = 1:InPatch_Residency{k}.EncounterNum
            range = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
            enc_bout_pix = [enc_bout_pix; round([nosetip(range,1) nosetip(range,2)])];
        end
        p_bout_img = (zeros(2040,2050));
        
        
        p_pix_indx = enc_bout_pix(:,2)+((enc_bout_pix(:,1)-1)*size(p_bout_img,1)); % linear index into bw
        p_pix_indx = p_pix_indx(isfinite(p_pix_indx));
        p_bout_img(p_pix_indx) = 1;
        
        ttt = bwmorph(p_bout_img,'bridge');
        ttt = bwmorph(ttt,'diag');
        ttt = bwmorph(ttt,'dilate');
        ttt = bwmorph(ttt,'dilate');
        ttt = bwmorph(ttt,'dilate');
        ttt = bwmorph(ttt,'erode');
        ttt = bwmorph(ttt,'erode');
        ttt = bwmorph(ttt,'majority');
        ttt = bwmorph(ttt,'close');
        
        [L,NUM] = bwlabel(ttt);
        STATS = regionprops(L,{'Area','PixelList'});
        indx_2 = find([STATS.Area] == max([STATS.Area]));indx_2 = indx_2(1);
        
        InPatch_Residency{k}.TotalEncounterArea = STATS(indx_2).PixelList;
        match_pix = ismember(tmp_patch_pixels,STATS(indx_2).PixelList,'rows');
        
        InPatch_Residency{k}.TotalEncounterAreaCovered = (sum(match_pix)/length(match_pix))*100;
    else
        InPatch_Residency{k}.TotalEncounterAreaCovered = 0;
    end
    %%%%%%%%%%%%%%%%
    
    feed_frames = {};
    feed_1stframes = [];
    for j = 1:InPatch_Residency{k}.Num
        feed_frames{j} = [InPatch_Residency{k}.Frames(j,1):InPatch_Residency{k}.Frames(j,2)];
        feed_1stframes(j) = InPatch_Residency{k}.Frames(j,1);
    end
    
    for i = 1:InPatch_Residency{k}.EncounterNum
        flag = 0;
        for j = 1:InPatch_Residency{k}.Num
            if ismember(InPatch_Residency{k}.EncounterFrames(i,1),feed_frames{j})
                InPatch_Residency{k}.Encounter_FeedingIndx(i) = 1;
                InPatch_Residency{k}.Encounter_FeedingNum(i) = j;
                InPatch_Residency{k}.Encounter_FoodPerc(i) = 100;
                flag = 1;
            end
        end
        if ~flag
            InPatch_Residency{k}.Encounter_FeedingIndx(i) = 0;
            tmp = (feed_1stframes-InPatch_Residency{k}.EncounterFrames(i,1));
            tmp_indx = tmp < 0;
            [tmp_find tmp_find_indx] = (max(tmp(tmp_indx)));
           
            
            if ~isempty(tmp_find_indx)
                InPatch_Residency{k}.Encounter_FeedingNum(i) = (tmp_find_indx);
                
                curr_feed_area = [];
                for q = 1:tmp_find_indx
                    curr_feed_area = [curr_feed_area; InPatch_Residency{k}.FeedingAreas{q}];
                end
                
                range = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
                encounter_pix = round([nosetip(range,1) nosetip(range,2)]);

                match_pix = ismember(encounter_pix,curr_feed_area,'rows');
                
                InPatch_Residency{k}.Encounter_FoodPerc(i) = 100-((sum(match_pix)/length(match_pix))*100);
            else
                InPatch_Residency{k}.Encounter_FeedingNum(i) = 0;
                InPatch_Residency{k}.Encounter_FoodPerc(i) = 0;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;hold on;
% for k = 1:NumPatches
%     
%     if InPatch_Residency{k}.EncounterNum ~= 0
%     plot(InPatch_Residency{k}.TotalEncounterArea(:,1),InPatch_Residency{k}.TotalEncounterArea(:,2),'.','color',[0.5 0.5 0.5]);
%     end
%     
%     C2 = colormap(jet(InPatch_Residency{k}.Num));
%     for i = 1:InPatch_Residency{k}.Num
%         plot(InPatch_Residency{k}.FeedingAreas{i}(:,1),InPatch_Residency{k}.FeedingAreas{i}(:,2),'.','color',C2(i,:));
%     end
%     a = PatchEdges{k};
%     plot(a(:,1),a(:,2),'k.')
%     
% %     C2 = colormap(jet(InPatch_Residency{k}.EncounterNum));
% %     for i = 1:InPatch_Residency{k}.EncounterNum
% %         if InPatch_Residency{k}.Encounter_FeedingIndx(i) == 0
% %         range = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
% %         plot(nosetip(range,1),nosetip(range,2),'.','color',C2(i,:));
% %         end
% %     end
% end
% axis equal
% 
% figure;hold on;
% for k = 1:NumPatches
%     C2 = colormap(jet(InPatch_Residency{k}.Num));
%     for i = 1:(InPatch_Residency{k}.Num)
%         range = InPatch_Residency{k}.Frames(i,1):InPatch_Residency{k}.Frames(i,2);
%         plot(nosetip(range,1),nosetip(range,2),'.','color',C2(i,:));
% 
% %         for j = InPatch_Residency{k}.Frames(i,1):InPatch_Residency{k}.Frames(i,2)
% %             plot(nosetailRT{j}(1,1),nosetailRT{j}(1,2),'.','color',C2(i,:));
% %         end
%     end
% end
% axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%
TotalFeedingTime = []; c10 = 0;

PercOfPatchEaten = [];
PercOfPatchVisited = [];

NumberOfFeedingsPerPatch = [];
PercentOfReturnsFreshEncounters = []; c9 = 0;
PercentOfReturnsFreshEncounters_Cumulative = [];

FeedingTimePerVisit = []; c3 = 0;
PatchCoveragePerVisit = [];

FeedingTimePerPatchVisit = []; c11 = 0;
PatchCoveragePerPatchVisit = [];

FeedingSearch_Feeding_Times = []; c6 = 0;
FeedingBreak_Feeding_Times = []; c7 = 0;

PercentOfPatchEatenWhenLeaving = []; c12 = 0;

NewFeedingVisitTime = []; c4 = 0;
NewFeedingVisitCoverage = [];

FirstFeedingVisitTime = []; 
FirstFeedingVisitCoverage = [];

ReturnFeedingVisitTime = []; c5 = 0;
ReturnFeedingVisitCoverage = [];

NonFeedingPatchEncounter_FoodEaten = []; c1 = 0;
NonFeedingPatchEncounter_CoM = [];
FeedingPatchEncounter_CoM = []; c2 = 0;

for k = 2:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(k) == 1
        c12 = c12+1;
        PatchLeft = PatchBouts.InPatchID(k-1);
        PatchLeftIndx = PatchBouts.InPatchIndx(k-1,1);
        PercentOfPatchEatenWhenLeaving(c12) = InPatch_Residency{PatchLeft}.CumPercEaten(find(InPatch_Residency{PatchLeft}.Frames(:,1)==PatchLeftIndx));
    end
end

for k = 1:size(InPatch_Residency,1)
    
    PercOfPatchEaten(k) = InPatch_Residency{k}.TotPercEaten;
    PercOfPatchVisited(k) = InPatch_Residency{k}.TotalEncounterAreaCovered;
    NumberOfFeedingsPerPatch(k) = InPatch_Residency{k}.Num;
    
    if InPatch_Residency{k}.Num > 1
        c9 = c9 + 1;
        PercentOfReturnsFreshEncounters(c9) = (sum(InPatch_Residency{k}.NewPatch(2:end))/(InPatch_Residency{k}.Num-1)) * 100;
        PercentOfReturnsFreshEncounters_Cumulative = [PercentOfReturnsFreshEncounters_Cumulative InPatch_Residency{k}.NewPatch(2:end)];
    end
    
    if InPatch_Residency{k}.Num > 0
        c10 = c10 + 1;
        TotalFeedingTime(c10) = sum(InPatch_Residency{k}.Durs);
        
        FirstFeedingVisitTime(c10) = InPatch_Residency{k}.Durs(1);
        FirstFeedingVisitCoverage(c10) = InPatch_Residency{k}.PercEaten(1);
        
        tmp_indx = [1:(InPatch_Residency{k}.Num)];
        tmp_indx2 = tmp_indx(logical(InPatch_Residency{k}.NewPatch));
        if length(tmp_indx2) == 1
            c11 = c11 + 1;
            temp1 = [];
            
            for ii = tmp_indx2(1):InPatch_Residency{k}.Num
                temp1 = [temp1 InPatch_Residency{k}.Durs(ii)];
            end
            FeedingTimePerPatchVisit(c11) = sum(temp1);
            PatchCoveragePerPatchVisit(c11) = InPatch_Residency{k}.TotPercEaten;
        else
            
            for i = 1:sum(InPatch_Residency{k}.NewPatch)-1
                c11 = c11 + 1;
                temp1 = [];
                
                for ii = tmp_indx2(i):tmp_indx2(i+1)-1
                    temp1 = [temp1 InPatch_Residency{k}.Durs(ii)];
                end
                FeedingTimePerPatchVisit(c11) = sum(temp1);
                if i == 1
                    PatchCoveragePerPatchVisit(c11) = InPatch_Residency{k}.CumPercEaten(tmp_indx2(i+1)-1);
                else
                    PatchCoveragePerPatchVisit(c11) = InPatch_Residency{k}.CumPercEaten(tmp_indx2(i+1)-1)-InPatch_Residency{k}.CumPercEaten(tmp_indx2(i-1));
                end
            end
            c11 = c11 + 1;
            temp1 = [];
            for ii = tmp_indx2(end):InPatch_Residency{k}.Num
                temp1 = [temp1 InPatch_Residency{k}.Durs(ii)];
            end
            FeedingTimePerPatchVisit(c11) = sum(temp1);
            PatchCoveragePerPatchVisit(c11) = InPatch_Residency{k}.CumPercEaten(end)-InPatch_Residency{k}.CumPercEaten(tmp_indx2(end)-1);
            
        end
    end
    
    TotalFeedingTime_Temp = [];
    for i = 1:InPatch_Residency{k}.Num
        c3 = c3 + 1;
        range = InPatch_Residency{k}.Frames(i,1):InPatch_Residency{k}.Frames(i,2);
        FeedingTimePerVisit(c3) = InPatch_Residency{k}.Durs(i);
        PatchCoveragePerVisit(c3) = InPatch_Residency{k}.PercEaten(i);
        
        if InPatch_Residency{k}.NewPatch(i) == 1
            c4 = c4 + 1;
            NewFeedingVisitTime(c4) = InPatch_Residency{k}.Durs(i);
            NewFeedingVisitCoverage(c4) = InPatch_Residency{k}.PercEaten(i);
            
            temp = InPatch_Residency{k}.Frames(i,1);
            temp2 = find(PatchBouts.InPatchIndx(:,1) == temp);
            temp2 = temp2 - 1;
            if temp2 ~= 0
                range1 = PatchBouts.InPatchIndx(temp2,2):InPatch_Residency{k}.Frames(i,1);
                
                c6 = c6 + 1;
                FeedingSearch_Feeding_Times(c6,1) = length(range1)/frameRate;
                FeedingSearch_Feeding_Times(c6,2) = InPatch_Residency{k}.Durs(i);
                FeedingSearch_Feeding_Times(c6,3) = InPatch_Residency{k}.PercEaten(i);
            end
            
        else
            c5 = c5 + 1;
            ReturnFeedingVisitTime(c5) = InPatch_Residency{k}.Durs(i);
            ReturnFeedingVisitCoverage(c5) = InPatch_Residency{k}.PercEaten(i);
            
            temp = InPatch_Residency{k}.Frames(i,1);
            temp2 = find(PatchBouts.InPatchIndx(:,1) == temp);
            temp2 = temp2 - 1;
            if temp2 ~= 0
                c7 = c7 + 1;
                range1 = PatchBouts.InPatchIndx(temp2,2):InPatch_Residency{k}.Frames(i,1);
                FeedingBreak_Feeding_Times(c7,1) = length(range1)/frameRate;
                FeedingBreak_Feeding_Times(c7,2) = InPatch_Residency{k}.Durs(i);
                FeedingBreak_Feeding_Times(c7,3) = InPatch_Residency{k}.PercEaten(i);
            end
        end
    end
    
    
    
    for i = 1:InPatch_Residency{k}.EncounterNum
        if InPatch_Residency{k}.Encounter_FeedingIndx(i) == 0
            c1 = c1 + 1;
            range = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
            NonFeedingPatchEncounter_CoM(c1) = nanmean(CoM_Vel(range));
            NonFeedingPatchEncounter_FoodEaten(c1) = InPatch_Residency{k}.Encounter_FoodPerc(i);
        else
            c2 = c2 + 1;
            range = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
            FeedingPatchEncounter_CoM(c2) = nanmean(CoM_Vel(range));
        end
    end
end

%%

save([Worm ' PatchAnalysis'],'Worm','nosetip','NumPatches','PatchResidency','NoseToPatchDistances','t_WL','Closeness',...
    'BinaryImage','PatchCentroids','PatchImages','PatchPixels','PatchEdges','PatchMajAxis','PatchAreas',...
    'CoM_Vel','CoM_Vel_Roam','PatchBouts','InPatch_Residency','TotalFeedingTime','PercOfPatchEaten',...
    'PercOfPatchVisited','NumberOfFeedingsPerPatch','PercentOfReturnsFreshEncounters','PercentOfReturnsFreshEncounters_Cumulative','FeedingTimePerVisit','FeedingTimePerPatchVisit',...
    'PatchCoveragePerVisit','PatchCoveragePerPatchVisit','FeedingSearch_Feeding_Times','FeedingBreak_Feeding_Times','NewFeedingVisitTime',...
    'NewFeedingVisitCoverage','FirstFeedingVisitTime','FirstFeedingVisitCoverage','ReturnFeedingVisitTime','ReturnFeedingVisitCoverage',...
    'NonFeedingPatchEncounter_FoodEaten','NonFeedingPatchEncounter_CoM','FeedingPatchEncounter_CoM');







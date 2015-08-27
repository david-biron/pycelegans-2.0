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

%%
y = AllHeadAnglesRT(1,:);
AnteriorAngles = zeros(18,length(y)-10);
AnteriorAngles_smooth = zeros(18,length(y)-10);
for j = 1:18
    y = AllHeadAnglesRT(j,:);
    tt = 0;
    ttt = 0;
    for i = 11:length(y)-10
        t = nanmean(y(i-10:i+10));
        if ~isfinite(t)
            t = tt;
        else
            tt = t;
        end
        if isfinite(y(i))
            AnteriorAngles(j,i) = y(i)-t;
            ttt = AnteriorAngles(j,i);
        else
            AnteriorAngles(j,i) =  ttt;
        end
    end
    AnteriorAngles_smooth(j,:) = smooth(AnteriorAngles(j,:),5);
end

AnteriorBendPeaks_Pos_x = {};
AnteriorBendPeaks_Pos_y = {};
AnteriorBendPeaks_Neg_x = {};
AnteriorBendPeaks_Neg_y = {};
Ant_Bend_Pos = zeros(18,size(RT_L,2))*NaN;
Ant_Bend_Neg = zeros(18,size(RT_L,2))*NaN;
for j = 1:18
    [AnteriorBendPeaks_Pos_x{j},AnteriorBendPeaks_Pos_y{j}] = peakfinder(AnteriorAngles_smooth(j,:),0.2);
    [AnteriorBendPeaks_Neg_x{j},AnteriorBendPeaks_Neg_y{j}] = peakfinder(-AnteriorAngles_smooth(j,:),0.2);
    
    
    Ant_Bend_Pos(j,AnteriorBendPeaks_Pos_x{j}) = AnteriorBendPeaks_Pos_y{j};
    Ant_Bend_Neg(j,AnteriorBendPeaks_Neg_x{j}) = AnteriorBendPeaks_Neg_y{j};
end

%%%%%%%%%%
AnteriorAngleRatio_1min = zeros(18,size(RT_L,2))*NaN;
AnteriorAngleAmplitude_1min = zeros(18,size(RT_L,2))*NaN;
for j = 1:18
    for i = 26:size(RT_L,2)-24
        t = nanmean(Ant_Bend_Pos(j,i-25:i+24))/nanmean(Ant_Bend_Neg(j,i-25:i+24))-1;
        t(t>1) = NaN;
        t(t<-1) = NaN;
        AnteriorAngleRatio_1min(j,i) = t;
        
        AnteriorAngleAmplitude_1min(j,i) = (nanmean(Ant_Bend_Pos(j,i-25:i+24))+nanmean(Ant_Bend_Neg(j,i-25:i+24)))/2;
    end
end

AnteriorAngleAmplitude_C_new = {};c1 = 0;
AnteriorAngleAmplitude_C_return = {};c2 = 0;
AnteriorAngleAmplitude_C_leave = {};c3 = 0;

AnteriorAngleRatio_C_new = {};
AnteriorAngleRatio_C_return = {};
AnteriorAngleRatio_C_leave = {};

for i = 1:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            c1 = c1 + 1;
            for j = 1:18
                range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
                if length(range1) > 3000
                    AnteriorAngleAmplitude_C_new{j}(c1,1:3000) = abs(AnteriorAngleAmplitude_1min(j,range1(end-2999:end)));
                    AnteriorAngleRatio_C_new{j}(c1,1:3000) = abs(AnteriorAngleRatio_1min(j,range1(end-2999:end)));
                else
                    AnteriorAngleAmplitude_C_new{j}(c1,1:3000-length(range1)) = NaN;
                    AnteriorAngleAmplitude_C_new{j}(c1,3000-length(range1)+1:3000) = abs(AnteriorAngleAmplitude_1min(j,range1));
                    
                    AnteriorAngleRatio_C_new{j}(c1,1:3000-length(range1)) = NaN;
                    AnteriorAngleRatio_C_new{j}(c1,3000-length(range1)+1:3000) = abs(AnteriorAngleRatio_1min(j,range1));
                end
                range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
                if length(range2) > 3000
                    AnteriorAngleAmplitude_C_new{j}(c1,3001:6000) = abs(AnteriorAngleAmplitude_1min(j,range2(1:3000)));
                    
                    AnteriorAngleRatio_C_new{j}(c1,3001:6000) = abs(AnteriorAngleRatio_1min(j,range2(1:3000)));
                else
                    AnteriorAngleAmplitude_C_new{j}(c1,3001:3000+length(range2)) = abs(AnteriorAngleAmplitude_1min(j,range2));
                    AnteriorAngleAmplitude_C_new{j}(c1,6000-(3000-length(range2)):6000) = NaN;
                    
                    AnteriorAngleRatio_C_new{j}(c1,3001:3000+length(range2)) = abs(AnteriorAngleRatio_1min(j,range2));
                    AnteriorAngleRatio_C_new{j}(c1,6000-(3000-length(range2)):6000) = NaN;
                end
            end
            
            c3 = c3 + 1;
            for j = 1:18
                range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
                if length(range1) > 3000
                    AnteriorAngleAmplitude_C_leave{j}(c1,1:3000) = abs(AnteriorAngleAmplitude_1min(j,range1(end-2999:end)));
                    
                    AnteriorAngleRatio_C_leave{j}(c1,1:3000) = abs(AnteriorAngleRatio_1min(j,range1(end-2999:end)));
                else
                    AnteriorAngleAmplitude_C_leave{j}(c1,1:3000-length(range1)) = NaN;
                    AnteriorAngleAmplitude_C_leave{j}(c1,3000-length(range1)+1:3000) = abs(AnteriorAngleAmplitude_1min(j,range1));
                    
                    AnteriorAngleRatio_C_leave{j}(c1,1:3000-length(range1)) = NaN;
                    AnteriorAngleRatio_C_leave{j}(c1,3000-length(range1)+1:3000) = abs(AnteriorAngleRatio_1min(j,range1));
                end
                range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
                if length(range2) > 3000
                    AnteriorAngleAmplitude_C_leave{j}(c1,3001:6000) = abs(AnteriorAngleAmplitude_1min(j,range2(1:3000)));
                    
                    AnteriorAngleRatio_C_leave{j}(c1,3001:6000) = abs(AnteriorAngleRatio_1min(j,range2(1:3000)));
                else
                    AnteriorAngleAmplitude_C_leave{j}(c1,3001:3000+length(range2)) = abs(AnteriorAngleAmplitude_1min(j,range2));
                    AnteriorAngleAmplitude_C_leave{j}(c1,6000-(3000-length(range2)):6000) = NaN;
                    
                    AnteriorAngleRatio_C_leave{j}(c1,3001:3000+length(range2)) = abs(AnteriorAngleRatio_1min(j,range2));
                    AnteriorAngleRatio_C_leave{j}(c1,6000-(3000-length(range2)):6000) = NaN;
                end
            end
        end
    else
        c2 = c2 + 1;
        for j = 1:18
            range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
            if length(range1) > 3000
                AnteriorAngleAmplitude_C_return{j}(c2,1:3000) = abs(AnteriorAngleAmplitude_1min(j,range1(end-2999:end)));
                
                AnteriorAngleRatio_C_return{j}(c2,1:3000) = abs(AnteriorAngleRatio_1min(j,range1(end-2999:end)));
            else
                AnteriorAngleAmplitude_C_return{j}(c2,1:3000-length(range1)) = NaN;
                AnteriorAngleAmplitude_C_return{j}(c2,3000-length(range1)+1:3000) = abs(AnteriorAngleAmplitude_1min(j,range1));
                
                AnteriorAngleRatio_C_return{j}(c2,1:3000-length(range1)) = NaN;
                AnteriorAngleRatio_C_return{j}(c2,3000-length(range1)+1:3000) = abs(AnteriorAngleRatio_1min(j,range1));
            end
            range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
            if length(range2) > 3000
                AnteriorAngleAmplitude_C_return{j}(c2,3001:6000) = abs(AnteriorAngleAmplitude_1min(j,range2(1:3000)));
                
                AnteriorAngleRatio_C_return{j}(c2,3001:6000) = abs(AnteriorAngleRatio_1min(j,range2(1:3000)));
            else
                AnteriorAngleAmplitude_C_return{j}(c2,3001:3000+length(range2)) = abs(AnteriorAngleAmplitude_1min(j,range2));
                AnteriorAngleAmplitude_C_return{j}(c2,6000-(3000-length(range2)):6000) = NaN;
                
                AnteriorAngleRatio_C_return{j}(c2,3001:3000+length(range2)) = abs(AnteriorAngleRatio_1min(j,range2));
                AnteriorAngleRatio_C_return{j}(c2,6000-(3000-length(range2)):6000) = NaN;
            end
        end
    end
end

AnteriorAngleAmplitude_new_Avg=[];AnteriorAngleAmplitude_new_STD=[];AnteriorAngleAmplitude_new_Ns=[];
AnteriorAngleAmplitude_return_Avg=[];AnteriorAngleAmplitude_return_STD=[];AnteriorAngleAmplitude_return_Ns=[];
AnteriorAngleAmplitude_leave_Avg=[];AnteriorAngleAmplitude_leave_STD=[];AnteriorAngleAmplitude_leave_Ns=[];

for j = 1:18
    [AnteriorAngleAmplitude_new_Avg(j,:),AnteriorAngleAmplitude_new_STD(j,:),AnteriorAngleAmplitude_new_Ns(j,:)] = av_fins(AnteriorAngleAmplitude_C_new{j});
    [AnteriorAngleAmplitude_return_Avg(j,:),AnteriorAngleAmplitude_return_STD(j,:),AnteriorAngleAmplitude_return_Ns(j,:)] = av_fins(AnteriorAngleAmplitude_C_return{j});
    [AnteriorAngleAmplitude_leave_Avg(j,:),AnteriorAngleAmplitude_leave_STD(j,:),AnteriorAngleAmplitude_leave_Ns(j,:)] = av_fins(AnteriorAngleAmplitude_C_leave{j});
end

AnteriorAngleRatio_new_Avg=[];AnteriorAngleRatio_new_STD=[];AnteriorAngleRatio_new_Ns=[];
AnteriorAngleRatio_return_Avg=[];AnteriorAngleRatio_return_STD=[];AnteriorAngleRatio_return_Ns=[];
AnteriorAngleRatio_leave_Avg=[];AnteriorAngleRatio_leave_STD=[];AnteriorAngleRatio_leave_Ns=[];

for j = 1:18
    [AnteriorAngleRatio_new_Avg(j,:),AnteriorAngleRatio_new_STD(j,:),AnteriorAngleRatio_new_Ns(j,:)] = av_fins(AnteriorAngleRatio_C_new{j});
    [AnteriorAngleRatio_return_Avg(j,:),AnteriorAngleRatio_return_STD(j,:),AnteriorAngleRatio_return_Ns(j,:)] = av_fins(AnteriorAngleRatio_C_return{j});
    [AnteriorAngleRatio_leave_Avg(j,:),AnteriorAngleRatio_leave_STD(j,:),AnteriorAngleRatio_leave_Ns(j,:)] = av_fins(AnteriorAngleRatio_C_leave{j});
end
% 
% 
% x = [1:6000];x = x/frameRate;
% 
% C = colormap(jet(18));
% figure;hold on;
% for j = 1:18
%     plot(x,AnteriorAngleAmplitude_new_Avg(j,:),'linewidth',2,'color',C(j,:))
%     % plot(x,AnteriorAngleAmplitude_new_Avg+AnteriorAngleAmplitude_new_STD./sqrt(AnteriorAngleAmplitude_new_Ns),'linewidth',1,'color',clr1*0.9)
%     % plot(x,AnteriorAngleAmplitude_new_Avg-AnteriorAngleAmplitude_new_STD./sqrt(AnteriorAngleAmplitude_new_Ns),'linewidth',1,'color',clr1*1.1)
% end
% 
% figure;hold on;
% for j = 1:18
% plot(x,AnteriorAngleAmplitude_return_Avg(j,:),'linewidth',2,'color',C(j,:))
% % plot(x,AnteriorAngleAmplitude_return_Avg+AnteriorAngleAmplitude_return_STD./sqrt(AnteriorAngleAmplitude_return_Ns),'linewidth',1,'color',clr1*0.9)
% % plot(x,AnteriorAngleAmplitude_return_Avg-AnteriorAngleAmplitude_return_STD./sqrt(AnteriorAngleAmplitude_return_Ns),'linewidth',1,'color',clr1*1.1)
% end
% 
% figure;hold on;
% for j = 1:18
% plot(x,AnteriorAngleAmplitude_leave_Avg(j,:),'linewidth',2,'color',C(j,:))
% % plot(x,AnteriorAngleAmplitude_leave_Avg+AnteriorAngleAmplitude_leave_STD./sqrt(AnteriorAngleAmplitude_leave_Ns),'linewidth',1,'color',clr1*0.9)
% % plot(x,AnteriorAngleAmplitude_leave_Avg-AnteriorAngleAmplitude_leave_STD./sqrt(AnteriorAngleAmplitude_leave_Ns),'linewidth',1,'color',clr1*1.1)
% end

% %%%%%
% 
% C = colormap(jet(18));
% figure;hold on;
% for j = 1:18
%     plot(x,AnteriorAngleRatio_new_Avg(j,:),'linewidth',1,'color',C(j,:))
% end
% 
% figure;hold on;
% for j = 1:18
% plot(x,AnteriorAngleRatio_return_Avg(j,:),'linewidth',1,'color',C(j,:))
% end
% 
% figure;hold on;
% for j = 1:18
% plot(x,AnteriorAngleRatio_leave_Avg(j,:),'linewidth',1,'color',C(j,:))
% end
%%
PrimaryBehavior = nan(1,length(RT_L));
for k = 1:length(RT_L)
    if isfinite(Behavior(:,k))
        PrimaryBehavior(k) = find(Behavior(:,k) == max(Behavior(:,k)),1,'first');
    end
end
    
indx = isfinite(PrimaryBehavior);
PrimaryBehaviorTemp = PrimaryBehavior(indx);
    
timeTemp = [1:length(PrimaryBehavior)];
timeTemp = timeTemp(indx);
    
BehaviorEdges_t = find(diff(PrimaryBehaviorTemp));
BehaviorEdges = timeTemp(BehaviorEdges_t);
%%%%%%%%%%%%%%%%%%%
BehaviorEdges_Logical = zeros(1,size(RT_L,2));
BehaviorEdges_Logical(BehaviorEdges) = 1;

BehaviorSwitches = zeros(1,size(RT_L,2))*NaN;

for i = 76:size(RT_L,2)-74
    BehaviorSwitches(i) = sum(BehaviorEdges_Logical(i-75:i+74));
end

BehSwitches_C_new = [];c1 = 0;
BehSwitches_C_return = [];c2 = 0;
BehSwitches_C_leave = [];c3 = 0;
for i = 1:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            c1 = c1 + 1;
            range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
            if length(range1) > 3000
                BehSwitches_C_new(c1,1:3000) = BehaviorSwitches(range1(end-2999:end));
            else
                BehSwitches_C_new(c1,1:3000-length(range1)) = NaN;
                BehSwitches_C_new(c1,3000-length(range1)+1:3000) = BehaviorSwitches(range1);
            end
            range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
            if length(range2) > 3000
                BehSwitches_C_new(c1,3001:6000) = BehaviorSwitches(range2(1:3000));
            else
                BehSwitches_C_new(c1,3001:3000+length(range2)) = BehaviorSwitches(range2);
                BehSwitches_C_new(c1,6000-(3000-length(range2)):6000) = NaN;
            end
            
            c3 = c3 + 1;
            range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
            if length(range1) > 3000
                BehSwitches_C_leave(c1,1:3000) = BehaviorSwitches(range1(end-2999:end));
            else
                BehSwitches_C_leave(c1,1:3000-length(range1)) = NaN;
                BehSwitches_C_leave(c1,3000-length(range1)+1:3000) = BehaviorSwitches(range1);
            end
            range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
            if length(range2) > 3000
                BehSwitches_C_leave(c1,3001:6000) = BehaviorSwitches(range2(1:3000));
            else
                BehSwitches_C_leave(c1,3001:3000+length(range2)) = BehaviorSwitches(range2);
                BehSwitches_C_leave(c1,6000-(3000-length(range2)):6000) = NaN;
            end
        end
    else
        c2 = c2 + 1;
        range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
        if length(range1) > 3000
            BehSwitches_C_return(c2,1:3000) = BehaviorSwitches(range1(end-2999:end));
        else
            BehSwitches_C_return(c2,1:3000-length(range1)) = NaN;
            BehSwitches_C_return(c2,3000-length(range1)+1:3000) = BehaviorSwitches(range1);
        end
        range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
        if length(range2) > 3000
            BehSwitches_C_return(c2,3001:6000) = BehaviorSwitches(range2(1:3000));
        else
            BehSwitches_C_return(c2,3001:3000+length(range2)) = BehaviorSwitches(range2);
            BehSwitches_C_return(c2,6000-(3000-length(range2)):6000) = NaN;
        end
    end
end
% 
[BehSwitches_new_Avg,BehSwitches_new_STD,BehSwitches_new_Ns] = av_fins(BehSwitches_C_new);
[BehSwitches_return_Avg,BehSwitches_return_STD,BehSwitches_return_Ns] = av_fins(BehSwitches_C_return);
[BehSwitches_leave_Avg,BehSwitches_leave_STD,BehSwitches_leave_Ns] = av_fins(BehSwitches_C_leave);
% 
% 
% x = [1:6000];x = x/frameRate;
% figure;hold on;
% plot(x,BehSwitches_new_Avg,'linewidth',2,'color',clr1)
% plot(x,BehSwitches_new_Avg+BehSwitches_new_STD./sqrt(BehSwitches_new_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,BehSwitches_new_Avg-BehSwitches_new_STD./sqrt(BehSwitches_new_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,BehSwitches_return_Avg,'linewidth',2,'color',clr1)
% plot(x,BehSwitches_return_Avg+BehSwitches_return_STD./sqrt(BehSwitches_return_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,BehSwitches_return_Avg-BehSwitches_return_STD./sqrt(BehSwitches_return_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,BehSwitches_leave_Avg,'linewidth',2,'color',clr1)
% plot(x,BehSwitches_leave_Avg+BehSwitches_leave_STD./sqrt(BehSwitches_leave_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,BehSwitches_leave_Avg-BehSwitches_leave_STD./sqrt(BehSwitches_leave_Ns),'linewidth',1,'color',clr1*1.1)
%%


NoseToPatchCentroids = zeros(NumPatches,size(nosetip,1))*NaN;
NoseToPatchCentroids_theta = zeros(NumPatches,size(nosetip,1))*NaN;

for j = 1:NumPatches
    t_dist = pdist2(nosetip(:,:),PatchCentroids(j,:));
    NoseToPatchCentroids(j,:) = t_dist;
    
    x = nosetip(:,1)-PatchCentroids(j,1);
    y = nosetip(:,2)-PatchCentroids(j,2);
    NoseToPatchCentroids_theta(j,:) = atan2(y , x);
end

DistFromCentroid_C_feeding = [];
DistFromCentroid_C_leave = []; DistFromCentroid_C_leave_NumEvents = 0;
ThetaFromCentroid_leave = [];

CoordsForLeavingEvents = {};
DistFromCentroidEvents = {};
ThetaForLeavingEvents = {};
LeavingEventsRadialRange = [];

for i = 1:size(PatchBouts.InPatchIndx,1)
    p_id = PatchBouts.InPatchID(i);
    range = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
    DistFromCentroid_C_feeding = [DistFromCentroid_C_feeding NoseToPatchCentroids(p_id,range)];
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            p_id = PatchBouts.InPatchID(i-1);
            
            range = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
            tmp = NoseToPatchCentroids(p_id,range);
            tmp_indx = find(tmp>t_WL*6);
            if ~isempty(tmp_indx)
                tmp_indx = tmp_indx(1);
                DistFromCentroid_C_leave_NumEvents = DistFromCentroid_C_leave_NumEvents+1;
                DistFromCentroid_C_leave = [DistFromCentroid_C_leave tmp(1:tmp_indx)];
                DistFromCentroidEvents{DistFromCentroid_C_leave_NumEvents} = tmp(1:tmp_indx);
                tmp2 = NoseToPatchCentroids_theta(p_id,range);
                ThetaFromCentroid_leave = [ThetaFromCentroid_leave  tmp2(1:tmp_indx)];
%                 
                CoordsForLeavingEvents{DistFromCentroid_C_leave_NumEvents} = nosetip(range(1:tmp_indx),:);
                
                ThetaForLeavingEvents{DistFromCentroid_C_leave_NumEvents} = NoseToPatchCentroids_theta(p_id,range(1:tmp_indx),:);
                
                tmp4 = tmp(1:tmp_indx);
                tmp4_indx = tmp4>40;
                tmp3 = radtodeg(ThetaForLeavingEvents{DistFromCentroid_C_leave_NumEvents});
                tmp3 = tmp3(tmp4_indx);
                bins = [-180:10:180];
                [tmpy tmpx] = hist(tmp3,bins);
                
                LeavingEventsRadialRange(DistFromCentroid_C_leave_NumEvents) = (sum(tmpy~=0)/size(bins,2)) * 100;
            end
        end
    end
end
% 
% bins = [0:5:ceil(t_WL*6)];
% [y1 x1] = hist(DistFromCentroid_C_leave,bins); y1 = y1./sum(y1);
% 
% [tF1 tGOF1] = fit(x1',y1','gauss3');
% figure;hold on;
% plot(x1,y1,'k');
% plot(tF1);


%%
%%%%%%%%%%%%%%%%%%%%%
DistFromPatch_C_new = [];c1 = 0;
DistFromPatch_C_return = [];c2 = 0;
DistFromPatch_C_leave = [];c3 = 0;
DistFromPatch_C_break = [];c4 = 0;

for i = 1:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            c1 = c1 + 1;
            p_id = PatchBouts.InPatchID(i);
            range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
            if length(range1) > 3000
                DistFromPatch_C_new(c1,1:3000) = NoseToPatchDistances(p_id,range1(end-2999:end));
            else
                DistFromPatch_C_new(c1,1:3000-length(range1)) = NaN;
                DistFromPatch_C_new(c1,3000-length(range1)+1:3000) = NoseToPatchDistances(p_id,range1);
            end
            range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
            if length(range2) > 3000
                DistFromPatch_C_new(c1,3001:6000) = NoseToPatchDistances(p_id,range2(1:3000));
            else
                DistFromPatch_C_new(c1,3001:3000+length(range2)) = NoseToPatchDistances(p_id,range2);
                DistFromPatch_C_new(c1,6000-(3000-length(range2)):6000) = NaN;
            end
            
            c3 = c3 + 1;
            p_id = PatchBouts.InPatchID(i-1);
            range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
            if length(range1) > 3000
                DistFromPatch_C_leave(c3,1:3000) = NoseToPatchDistances(p_id,range1(end-2999:end));
            else
                DistFromPatch_C_leave(c3,1:3000-length(range1)) = NaN;
                DistFromPatch_C_leave(c3,3000-length(range1)+1:3000) = NoseToPatchDistances(p_id,range1);
            end
            range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
            if length(range2) > 3000
                DistFromPatch_C_leave(c3,3001:6000) = NoseToPatchDistances(p_id,range2(1:3000));
            else
                DistFromPatch_C_leave(c3,3001:3000+length(range2)) = NoseToPatchDistances(p_id,range2);
                DistFromPatch_C_leave(c3,6000-(3000-length(range2)):6000) = NaN;
            end
        end
    else
        c2 = c2 + 1;
        p_id = PatchBouts.InPatchID(i);
        range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
        if length(range1) > 3000
            DistFromPatch_C_return(c2,1:3000) = NoseToPatchDistances(p_id,range1(end-2999:end));
        else
            DistFromPatch_C_return(c2,1:3000-length(range1)) = NaN;
            DistFromPatch_C_return(c2,3000-length(range1)+1:3000) = NoseToPatchDistances(p_id,range1);
        end
        range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
        if length(range2) > 3000
            DistFromPatch_C_return(c2,3001:6000) = NoseToPatchDistances(p_id,range2(1:3000));
        else
            DistFromPatch_C_return(c2,3001:3000+length(range2)) = NoseToPatchDistances(p_id,range2);
            DistFromPatch_C_return(c2,6000-(3000-length(range2)):6000) = NaN;
        end
        
        c4 = c4 + 1;
        p_id = PatchBouts.InPatchID(i-1);
        range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
        if length(range1) > 3000
            DistFromPatch_C_break(c4,1:3000) = NoseToPatchDistances(p_id,range1(end-2999:end));
        else
            DistFromPatch_C_break(c4,1:3000-length(range1)) = NaN;
            DistFromPatch_C_break(c4,3000-length(range1)+1:3000) = NoseToPatchDistances(p_id,range1);
        end
        range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
        if length(range2) > 3000
            DistFromPatch_C_break(c4,3001:6000) = NoseToPatchDistances(p_id,range2(1:3000));
        else
            DistFromPatch_C_break(c4,3001:3000+length(range2)) = NoseToPatchDistances(p_id,range2);
            DistFromPatch_C_break(c4,6000-(3000-length(range2)):6000) = NaN;
        end
    end
end

TimeToGetAway = [];
TimeToGetAway_long = []; c1 = 0;
TimeToGetAway_Coll = [];
for i = 1:size(DistFromPatch_C_leave,1)
    for j = 1:50
        tmp_y = DistFromPatch_C_leave(i,3000:end);
        tmp = find(tmp_y > t_WL*j/10);
        if ~isempty(tmp)
            tmp = tmp(1); tmp = (tmp)/frameRate;
        else
            tmp = NaN;
        end
        TimeToGetAway(j,i) = tmp;
    end
    tmp_y = DistFromPatch_C_leave(i,3000:end);
    if max(tmp_y(isfinite(tmp_y))) > t_WL * 3
        c1 = c1 + 1;
        TimeToGetAway_Coll = [TimeToGetAway_Coll tmp_y];
        for j = 1:50
            tmp = find(tmp_y > t_WL*j/10);
            if ~isempty(tmp)
                tmp = tmp(1); tmp = (tmp)/frameRate;
            else
                tmp = NaN;
            end
            TimeToGetAway_long(j,c1) = tmp;
        end
    end
end
TimeToReachDistanceFromPatch = [];
TimeToReachDistanceFromPatch(:,1) = [1:50]*10;
TimeToReachDistanceFromPatch(:,2) = nanmean(TimeToGetAway');

% nknots = 2;
% sigma = 0;
% options = struct('KnotRemoval','none','sigma',sigma);
% pp = BSFK(TimeToReachDistanceFromPatch(:,1),TimeToReachDistanceFromPatch(:,2),2,nknots,[],options);
% yfit = ppval(pp,TimeToReachDistanceFromPatch(:,1));

% pp = BSFK([1:50]*10,TimeToGetAway(:,1),2,nknots,[],options);
% yfit = ppval(pp,[1:50]*10);
% 
% x = pp.breaks(2);
% y = ppval(pp,x);
% LocalSearchToRoamCrossOver = [x y];
% LocalSearchToRoamSlopes = [pp.coefs(:,1)];

% figure;
% plot([1:50]*10,TimeToGetAway(:,1),'k')
% hold on;
% plot([1:50]*10,yfit,'r');

% figure;
% plot(TimeToReachDistanceFromPatch(:,1),TimeToReachDistanceFromPatch(:,2),'k')
% hold on;
% plot(TimeToReachDistanceFromPatch(:,1),yfit,'r');
% set(gcf,'position',scrsz);
% saveas(gcf,[Worm  ' LocalSearchToRoamCrossOver.pdf']);
% close

TimeToGetAway_Break = [];
for i = 1:size(DistFromPatch_C_break,1)
    for j = 1:50
        tmp_y = DistFromPatch_C_break(i,3000:end);
        tmp = find(tmp_y > t_WL*j/10);
        if ~isempty(tmp)
            tmp = tmp(1); tmp = (tmp)/frameRate;
        else
            tmp = NaN;
        end
        TimeToGetAway_Break(j,i) = tmp;
    end
end

MaxDistanceForFeedingBreaks = [];
MeanDistanceForFeedingBreaks = [];

for i = 1:size(DistFromPatch_C_break,1)
    tmp_y = DistFromPatch_C_break(i,3000:end);
    tmp = max(tmp_y);
    MaxDistanceForFeedingBreaks(i) = tmp;
    MeanDistanceForFeedingBreaks(i) = nanmean(tmp_y);
end


[DistFromPatch_new_Avg,DistFromPatch_new_STD,DistFromPatch_new_Ns] = av_fins(DistFromPatch_C_new);
[DistFromPatch_return_Avg,DistFromPatch_return_STD,DistFromPatch_return_Ns] = av_fins(DistFromPatch_C_return);
[DistFromPatch_leave_Avg,DistFromPatch_leave_STD,DistFromPatch_leave_Ns] = av_fins(DistFromPatch_C_leave);
[DistFromPatch_break_Avg,DistFromPatch_break_STD,DistFromPatch_break_Ns] = av_fins(DistFromPatch_C_break);

ApproachToPatch(:,1) = [1:600]/frameRate;
ApproachToPatch(:,2) = fliplr(DistFromPatch_new_Avg(2401:3000));

% [ApproachToPatchFit(1) ApproachToPatchFit(2) ApproachToPatchFit(3) ApproachToPatchFit(4)] = logfit(ApproachToPatch(:,1),ApproachToPatch(:,2),'loglog');close; %%% slope / intercept / Mean square error / coefficient of determination

% 
% x = [1:6000];x = x/frameRate;
% figure;hold on;
% plot(x(1:600),fliplr(DistFromPatch_new_Avg(2401:3000)),'linewidth',2,'color',clr1)
% 
% figure;hold on;
% plot(x,DistFromPatch_new_Avg,'linewidth',2,'color',clr1)
% plot(x,DistFromPatch_new_Avg+DistFromPatch_new_STD./sqrt(DistFromPatch_new_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,DistFromPatch_new_Avg-DistFromPatch_new_STD./sqrt(DistFromPatch_new_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,DistFromPatch_return_Avg,'linewidth',2,'color',clr1)
% plot(x,DistFromPatch_return_Avg+DistFromPatch_return_STD./sqrt(DistFromPatch_return_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,DistFromPatch_return_Avg-DistFromPatch_return_STD./sqrt(DistFromPatch_return_Ns),'linewidth',1,'color',clr1*1.1)
% % 
% x = [1:6000];x = x/frameRate;
% figure;hold on;
% plot(x,DistFromPatch_leave_Avg,'linewidth',2,'color',clr1)
% plot(x,DistFromPatch_leave_Avg+DistFromPatch_leave_STD./sqrt(DistFromPatch_leave_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,DistFromPatch_leave_Avg-DistFromPatch_leave_STD./sqrt(DistFromPatch_leave_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,DistFromPatch_break_Avg,'linewidth',2,'color',clr1)
% plot(x,DistFromPatch_break_Avg+DistFromPatch_break_STD./sqrt(DistFromPatch_break_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,DistFromPatch_break_Avg-DistFromPatch_break_STD./sqrt(DistFromPatch_break_Ns),'linewidth',1,'color',clr1*1.1)

%%
CoM_V_C_new = [];c1 = 0;
CoM_V_C_return = [];c2 = 0;
CoM_V_C_leave = [];c3 = 0;
CoM_V_C_break = [];c4 = 0;

for i = 1:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            c1 = c1 + 1;
            range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
            if length(range1) > 300
                CoM_V_C_new(c1,1:300) = CoM_Vel(range1(end-299:end));
            else
                CoM_V_C_new(c1,1:300-length(range1)) = NaN;
                CoM_V_C_new(c1,300-length(range1)+1:300) = CoM_Vel(range1);
            end
            range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
            if length(range2) > 300
                CoM_V_C_new(c1,301:600) = CoM_Vel(range2(1:300));
            else
                CoM_V_C_new(c1,301:300+length(range2)) = CoM_Vel(range2);
                CoM_V_C_new(c1,600-(300-length(range2)):600) = NaN;
            end
            
            c3 = c3 + 1;
            range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
            if length(range1) > 300
                CoM_V_C_leave(c3,1:300) = CoM_Vel(range1(end-299:end));
            else
                CoM_V_C_leave(c3,1:300-length(range1)) = NaN;
                CoM_V_C_leave(c3,300-length(range1)+1:300) = CoM_Vel(range1);
            end
            range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
            if length(range2) > 300
                CoM_V_C_leave(c3,301:600) = CoM_Vel(range2(1:300));
            else
                CoM_V_C_leave(c3,301:300+length(range2)) = CoM_Vel(range2);
                CoM_V_C_leave(c3,600-(300-length(range2)):600) = NaN;
            end
        end
    else
        c2 = c2 + 1;
        range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
        if length(range1) > 300
            CoM_V_C_return(c2,1:300) = CoM_Vel(range1(end-299:end));
        else
            CoM_V_C_return(c2,1:300-length(range1)) = NaN;
            CoM_V_C_return(c2,300-length(range1)+1:300) = CoM_Vel(range1);
        end
        range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
        if length(range2) > 300
            CoM_V_C_return(c2,301:600) = CoM_Vel(range2(1:300));
        else
            CoM_V_C_return(c2,301:300+length(range2)) = CoM_Vel(range2);
            CoM_V_C_return(c2,600-(300-length(range2)):600) = NaN;
        end
        
        c4 = c4 + 1;
        range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
        if length(range1) > 300
            CoM_V_C_break(c4,1:300) = CoM_Vel(range1(end-299:end));
        else
            CoM_V_C_break(c4,1:300-length(range1)) = NaN;
            CoM_V_C_break(c4,300-length(range1)+1:300) = CoM_Vel(range1);
        end
        range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
        if length(range2) > 300
            CoM_V_C_break(c4,301:600) = CoM_Vel(range2(1:300));
        else
            CoM_V_C_break(c4,301:300+length(range2)) = CoM_Vel(range2);
            CoM_V_C_break(c4,600-(300-length(range2)):600) = NaN;
        end
    end
end

% [CoM_V_new_Avg,CoM_V_new_STD,CoM_V_new_Ns] = av_fins(CoM_V_C_new);
% [CoM_V_return_Avg,CoM_V_return_STD,CoM_V_return_Ns] = av_fins(CoM_V_C_return);
% [CoM_V_leave_Avg,CoM_V_leave_STD,CoM_V_leave_Ns] = av_fins(CoM_V_C_leave);
% [CoM_V_break_Avg,CoM_V_break_STD,CoM_V_break_Ns] = av_fins(CoM_V_C_break);
% 
% 
% x = [1:600];x = x/frameRate;
% figure;hold on;
% plot(x,CoM_V_new_Avg,'linewidth',2,'color',clr1)
% plot(x,CoM_V_new_Avg+CoM_V_new_STD./sqrt(CoM_V_new_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,CoM_V_new_Avg-CoM_V_new_STD./sqrt(CoM_V_new_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,CoM_V_return_Avg,'linewidth',2,'color',clr1)
% plot(x,CoM_V_return_Avg+CoM_V_return_STD./sqrt(CoM_V_return_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,CoM_V_return_Avg-CoM_V_return_STD./sqrt(CoM_V_return_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,CoM_V_leave_Avg,'linewidth',2,'color',clr1)
% plot(x,CoM_V_leave_Avg+CoM_V_leave_STD./sqrt(CoM_V_leave_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,CoM_V_leave_Avg-CoM_V_leave_STD./sqrt(CoM_V_leave_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,CoM_V_break_Avg,'linewidth',2,'color',clr1)
% plot(x,CoM_V_break_Avg+CoM_V_break_STD./sqrt(CoM_V_break_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,CoM_V_break_Avg-CoM_V_break_STD./sqrt(CoM_V_break_Ns),'linewidth',1,'color',clr1*1.1)

%%

CoM_Encounter_Feeding = [];c2 = 0;
CoM_Encounter_NonFeeding = [];c1 = 0;
CoM_Encounter_WithFood = [];c3 = 0;
CoM_Encounter_NoFood = [];c4 = 0;


for k = 1:size(InPatch_Residency,1)
    for i = 1:InPatch_Residency{k}.EncounterNum
        if InPatch_Residency{k}.Encounter_FeedingIndx(i) == 0
            temp = InPatch_Residency{k}.EncounterFrames(i,1);
            temp2 = find(PatchBouts.PatchCross(:,1) == temp);
            temp2 = temp2 - 1;
            if temp2 ~= 0
                c1 = c1 + 1;
                
                range1 = PatchBouts.PatchCross(temp2,2):(InPatch_Residency{k}.EncounterFrames(i,1)-1);
                if length(range1) > 300
                    CoM_Encounter_NonFeeding(c1,1:300) = CoM_Vel(range1(end-299:end));
                else
                    CoM_Encounter_NonFeeding(c1,1:300-length(range1)) = NaN;
                    CoM_Encounter_NonFeeding(c1,300-length(range1)+1:300) = CoM_Vel(range1);
                end
                range2 = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
                if length(range2) > 300
                    CoM_Encounter_NonFeeding(c1,301:600) = CoM_Vel(range2(1:300));
                else
                    CoM_Encounter_NonFeeding(c1,301:300+length(range2)) = CoM_Vel(range2);
                    CoM_Encounter_NonFeeding(c1,600-(300-length(range2)):600) = NaN;
                end
            end
            
            if InPatch_Residency{k}.Encounter_FoodPerc(i) < 30
                temp = InPatch_Residency{k}.EncounterFrames(i,1);
                temp2 = find(PatchBouts.PatchCross(:,1) == temp);
                temp2 = temp2 - 1;
                if temp2 ~= 0
                    c3 = c3 + 1;
                    
                    range1 = PatchBouts.PatchCross(temp2,2):(InPatch_Residency{k}.EncounterFrames(i,1)-1);
                    if length(range1) > 300
                        CoM_Encounter_WithFood(c3,1:300) = CoM_Vel(range1(end-299:end));
                    else
                        CoM_Encounter_WithFood(c3,1:300-length(range1)) = NaN;
                        CoM_Encounter_WithFood(c3,300-length(range1)+1:300) = CoM_Vel(range1);
                    end
                    range2 = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
                    if length(range2) > 300
                        CoM_Encounter_WithFood(c3,301:600) = CoM_Vel(range2(1:300));
                    else
                        CoM_Encounter_WithFood(c3,301:300+length(range2)) = CoM_Vel(range2);
                        CoM_Encounter_WithFood(c3,600-(300-length(range2)):600) = NaN;
                    end
                end
            end
            if InPatch_Residency{k}.Encounter_FoodPerc(i) > 70
                temp = InPatch_Residency{k}.EncounterFrames(i,1);
                temp2 = find(PatchBouts.PatchCross(:,1) == temp);
                temp2 = temp2 - 1;
                if temp2 ~= 0
                    c4 = c4 + 1;
                    
                    range1 = PatchBouts.PatchCross(temp2,2):(InPatch_Residency{k}.EncounterFrames(i,1)-1);
                    if length(range1) > 300
                        CoM_Encounter_NoFood(c4,1:300) = CoM_Vel(range1(end-299:end));
                    else
                        CoM_Encounter_NoFood(c4,1:300-length(range1)) = NaN;
                        CoM_Encounter_NoFood(c4,300-length(range1)+1:300) = CoM_Vel(range1);
                    end
                    range2 = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
                    if length(range2) > 300
                        CoM_Encounter_NoFood(c4,301:600) = CoM_Vel(range2(1:300));
                    else
                        CoM_Encounter_NoFood(c4,301:300+length(range2)) = CoM_Vel(range2);
                        CoM_Encounter_NoFood(c4,600-(300-length(range2)):600) = NaN;
                    end
                end
            end
        else
            
            temp = InPatch_Residency{k}.EncounterFrames(i,1);
            temp2 = find(PatchBouts.PatchCross(:,1) == temp);
            temp2 = temp2 - 1;
            if temp2 ~= 0
                c2 = c2 + 1;
                
                range1 = PatchBouts.PatchCross(temp2,2):(InPatch_Residency{k}.EncounterFrames(i,1)-1);
                if length(range1) > 300
                    CoM_Encounter_Feeding(c2,1:300) = CoM_Vel(range1(end-299:end));
                else
                    CoM_Encounter_Feeding(c2,1:300-length(range1)) = NaN;
                    CoM_Encounter_Feeding(c2,300-length(range1)+1:300) = CoM_Vel(range1);
                end
                range2 = InPatch_Residency{k}.EncounterFrames(i,1):InPatch_Residency{k}.EncounterFrames(i,2);
                if length(range2) > 300
                    CoM_Encounter_Feeding(c2,301:600) = CoM_Vel(range2(1:300));
                else
                    CoM_Encounter_Feeding(c2,301:300+length(range2)) = CoM_Vel(range2);
                    CoM_Encounter_Feeding(c2,600-(300-length(range2)):600) = NaN;
                end
                
            end
        end        
    end
end
% 
% [CoM_EncounterFeed_Avg,CoM_EncounterFeed_STD,CoM_EncounterFeed_Ns] = av_fins(CoM_Encounter_Feeding);
% [CoM_EncounterNonFeed_Avg,CoM_EncounterNonFeed_STD,CoM_EncounterNonFeed_Ns] = av_fins(CoM_Encounter_NonFeeding);
% [CoM_EncounterWithFood_Avg,CoM_EncounterWithFood_STD,CoM_EncounterWithFood_Ns] = av_fins(CoM_Encounter_WithFood);
% [CoM_EncounterWithoutFood_Avg,CoM_EncounterWithoutFood_STD,CoM_EncounterWithoutFood_Ns] = av_fins(CoM_Encounter_NoFood);
% 
% 
% x = [1:600];x = x/frameRate;
% figure;hold on;
% plot(x,CoM_EncounterFeed_Avg,'linewidth',2,'color',clr1)
% plot(x,CoM_EncounterFeed_Avg+CoM_EncounterFeed_STD./sqrt(CoM_EncounterFeed_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,CoM_EncounterFeed_Avg-CoM_EncounterFeed_STD./sqrt(CoM_EncounterFeed_Ns),'linewidth',1,'color',clr1*1.1)
% 
% % figure;hold on;
% plot(x,CoM_EncounterNonFeed_Avg,'linewidth',2,'color',clr2)
% plot(x,CoM_EncounterNonFeed_Avg+CoM_EncounterNonFeed_STD./sqrt(CoM_EncounterNonFeed_Ns),'linewidth',1,'color',clr2*0.9)
% plot(x,CoM_EncounterNonFeed_Avg-CoM_EncounterNonFeed_STD./sqrt(CoM_EncounterNonFeed_Ns),'linewidth',1,'color',clr2*1.1)
% 
% figure;hold on;
% plot(x,CoM_EncounterWithFood_Avg,'linewidth',2,'color',clr1)
% plot(x,CoM_EncounterWithFood_Avg+CoM_EncounterWithFood_STD./sqrt(CoM_EncounterWithFood_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,CoM_EncounterWithFood_Avg-CoM_EncounterWithFood_STD./sqrt(CoM_EncounterWithFood_Ns),'linewidth',1,'color',clr1*1.1)
% 
% % figure;hold on;
% plot(x,CoM_EncounterWithoutFood_Avg,'linewidth',2,'color',clr2)
% plot(x,CoM_EncounterWithoutFood_Avg+CoM_EncounterWithoutFood_STD./sqrt(CoM_EncounterWithoutFood_Ns),'linewidth',1,'color',clr2*0.9)
% plot(x,CoM_EncounterWithoutFood_Avg-CoM_EncounterWithoutFood_STD./sqrt(CoM_EncounterWithoutFood_Ns),'linewidth',1,'color',clr2*1.1)
%%
SumAngles_C_new = [];c1 = 0;
SumAngles_C_return = [];c2 = 0;
SumAngles_C_leave = [];c3 = 0;
for i = 1:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            c1 = c1 + 1;
            range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
            if length(range1) > 3000
                SumAngles_C_new(c1,1:3000) = sum(abs(AllHeadAnglesRT(:,range1(end-2999:end))));
            else
                SumAngles_C_new(c1,1:3000-length(range1)) = NaN;
                SumAngles_C_new(c1,3000-length(range1)+1:3000) = sum(abs(AllHeadAnglesRT(:,range1)));
            end
            range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
            if length(range2) > 3000
                SumAngles_C_new(c1,3001:6000) = sum(abs(AllHeadAnglesRT(:,range2(1:3000))));
            else
                SumAngles_C_new(c1,3001:3000+length(range2)) = sum(abs(AllHeadAnglesRT(:,range2)));
                SumAngles_C_new(c1,6000-(3000-length(range2)):6000) = NaN;
            end
            
            c3 = c3 + 1;
            range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
            if length(range1) > 3000
                SumAngles_C_leave(c1,1:3000) = sum(abs(AllHeadAnglesRT(:,range1(end-2999:end))));
            else
                SumAngles_C_leave(c1,1:3000-length(range1)) = NaN;
                SumAngles_C_leave(c1,3000-length(range1)+1:3000) = sum(abs(AllHeadAnglesRT(:,range1)));
            end
            range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
            if length(range2) > 3000
                SumAngles_C_leave(c1,3001:6000) = sum(abs(AllHeadAnglesRT(:,range2(1:3000))));
            else
                SumAngles_C_leave(c1,3001:3000+length(range2)) = sum(abs(AllHeadAnglesRT(:,range2)));
                SumAngles_C_leave(c1,6000-(3000-length(range2)):6000) = NaN;
            end
        end
    else
        c2 = c2 + 1;
        range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
        if length(range1) > 3000
            SumAngles_C_return(c2,1:3000) = sum(abs(AllHeadAnglesRT(:,range1(end-2999:end))));
        else
            SumAngles_C_return(c2,1:3000-length(range1)) = NaN;
            SumAngles_C_return(c2,3000-length(range1)+1:3000) = sum(abs(AllHeadAnglesRT(:,range1)));
        end
        range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
        if length(range2) > 3000
            SumAngles_C_return(c2,3001:6000) = sum(abs(AllHeadAnglesRT(:,range2(1:3000))));
        else
            SumAngles_C_return(c2,3001:3000+length(range2)) = sum(abs(AllHeadAnglesRT(:,range2)));
            SumAngles_C_return(c2,6000-(3000-length(range2)):6000) = NaN;
        end
    end
end
% 
% [SumAngles_new_Avg,SumAngles_new_STD,SumAngles_new_Ns] = av_fins(SumAngles_C_new);
% [SumAngles_return_Avg,SumAngles_return_STD,SumAngles_return_Ns] = av_fins(SumAngles_C_return);
% [SumAngles_leave_Avg,SumAngles_leave_STD,SumAngles_leave_Ns] = av_fins(SumAngles_C_leave);
% 
% 
% x = [1:6000];x = x/frameRate;
% figure;hold on;
% plot(x,SumAngles_new_Avg,'linewidth',2,'color',clr1)
% plot(x,SumAngles_new_Avg+SumAngles_new_STD./sqrt(SumAngles_new_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,SumAngles_new_Avg-SumAngles_new_STD./sqrt(SumAngles_new_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,SumAngles_return_Avg,'linewidth',2,'color',clr1)
% plot(x,SumAngles_return_Avg+SumAngles_return_STD./sqrt(SumAngles_return_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,SumAngles_return_Avg-SumAngles_return_STD./sqrt(SumAngles_return_Ns),'linewidth',1,'color',clr1*1.1)
% 
% figure;hold on;
% plot(x,SumAngles_leave_Avg,'linewidth',2,'color',clr1)
% plot(x,SumAngles_leave_Avg+SumAngles_leave_STD./sqrt(SumAngles_leave_Ns),'linewidth',1,'color',clr1*0.9)
% plot(x,SumAngles_leave_Avg-SumAngles_leave_STD./sqrt(SumAngles_leave_Ns),'linewidth',1,'color',clr1*1.1)
%%
dAngles_C_new = {};c1 = 0;
dAngles_C_return = {};c2 = 0;
dAngles_C_leave = {};c3 = 0;
for i = 1:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            c1 = c1 + 1;
            for j = 1:18
                range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
                if length(range1) > 300
                    dAngles_C_new{j}(c1,1:300) = abs(dAnglesRT(j,range1(end-299:end)));
                else
                    dAngles_C_new{j}(c1,1:300-length(range1)) = NaN;
                    dAngles_C_new{j}(c1,300-length(range1)+1:300) = abs(dAnglesRT(j,range1));
                end
                range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
                if length(range2) > 300
                    dAngles_C_new{j}(c1,301:600) = abs(dAnglesRT(j,range2(1:300)));
                else
                    dAngles_C_new{j}(c1,301:300+length(range2)) = abs(dAnglesRT(j,range2));
                    dAngles_C_new{j}(c1,600-(300-length(range2)):600) = NaN;
                end
            end
            for j = 1:18
                c3 = c3 + 1;
                range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
                if length(range1) > 300
                    dAngles_C_leave{j}(c1,1:300) = abs(dAnglesRT(j,range1(end-299:end)));
                else
                    dAngles_C_leave{j}(c1,1:300-length(range1)) = NaN;
                    dAngles_C_leave{j}(c1,300-length(range1)+1:300) = abs(dAnglesRT(j,range1));
                end
                range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
                if length(range2) > 300
                    dAngles_C_leave{j}(c1,301:600) = abs(dAnglesRT(j,range2(1:300)));
                else
                    dAngles_C_leave{j}(c1,301:300+length(range2)) = abs(dAnglesRT(j,range2));
                    dAngles_C_leave{j}(c1,600-(300-length(range2)):600) = NaN;
                end
            end
        end
    else
        c2 = c2 + 1;
        for j = 1:18
            range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
            if length(range1) > 300
                dAngles_C_return{j}(c2,1:300) = abs(dAnglesRT(j,range1(end-299:end)));
            else
                dAngles_C_return{j}(c2,1:300-length(range1)) = NaN;
                dAngles_C_return{j}(c2,300-length(range1)+1:300) = abs(dAnglesRT(j,range1));
            end
            range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
            if length(range2) > 300
                dAngles_C_return{j}(c2,301:600) = abs(dAnglesRT(j,range2(1:300)));
            else
                dAngles_C_return{j}(c2,301:300+length(range2)) = abs(dAnglesRT(j,range2));
                dAngles_C_return{j}(c2,600-(300-length(range2)):600) = NaN;
            end
        end
    end
end
% dAngles_new_Patch_Avg = []; dAngles_new_Patch_STD = []; dAngles_new_Patch_Ns = [];
% dAngles_return_Patch_Avg = []; dAngles_return_Patch_STD = []; dAngles_return_Patch_Ns = [];
% dAngles_leave_Patch_Avg = []; dAngles_leave_Patch_STD = []; dAngles_leave_Patch_Ns = [];
% 
% for j = 1:18
%     [dAngles_new_Patch_Avg(j,:), dAngles_new_Patch_STD(j,:), dAngles_new_Patch_Ns(j,:)] = av_fins(dAngles_C_new{j});
%     [dAngles_return_Patch_Avg(j,:), dAngles_return_Patch_STD(j,:), dAngles_return_Patch_Ns(j,:)] = av_fins(dAngles_C_return{j});
%     [dAngles_leave_Patch_Avg(j,:), dAngles_leave_Patch_STD(j,:), dAngles_leave_Patch_Ns(j,:)] = av_fins(dAngles_C_leave{j});
% end
% x = [1:600];x = x/frameRate;
% C = colormap(jet(18));
% figure;hold on;
% for j = 1:18
%    plot(x,dAngles_new_Patch_Avg(j,:),'color',C(j,:));
% end
% 
% figure;hold on;
% for j = 1:18
%    plot(x,dAngles_return_Patch_Avg(j,:),'color',C(j,:));
% end
% 
% figure;hold on;
% for j = 1:18
%    plot(x,dAngles_leave_Patch_Avg(j,:),'color',C(j,:));
% end

%%
%%%%%%%%%%%%%%%%%%
Behavior_C_new = {};c1 = 0;
Behavior_C_return = {};c2 = 0;
Behavior_C_leave = {};c3 = 0;
for i = 1:size(PatchBouts.InPatchIndx,1)
    if PatchBouts.NewPatchIndx(i) == 1
        if i == 1
        else
            c1 = c1 + 1;
            for j = 1:4
                range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
                if length(range1) > 300
                    Behavior_C_new{j}(c1,1:300) = Behavior(j,range1(end-299:end));
                else
                    Behavior_C_new{j}(c1,1:300-length(range1)) = NaN;
                    Behavior_C_new{j}(c1,300-length(range1)+1:300) = Behavior(j,range1);
                end
                range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
                if length(range2) > 300
                    Behavior_C_new{j}(c1,301:600) = Behavior(j,range2(1:300));
                else
                    Behavior_C_new{j}(c1,301:300+length(range2)) = Behavior(j,range2);
                    Behavior_C_new{j}(c1,600-(300-length(range2)):600) = NaN;
                end
            end
            
            c3 = c3 + 1;
            for j = 1:4
            range1 = PatchBouts.InPatchIndx(i-1,1):(PatchBouts.InPatchIndx(i-1,2)-1);
            if length(range1) > 300
                Behavior_C_leave{j}(c1,1:300) = Behavior(j,range1(end-299:end));
            else
                Behavior_C_leave{j}(c1,1:300-length(range1)) = NaN;
                Behavior_C_leave{j}(c1,300-length(range1)+1:300) = Behavior(j,range1);
            end
            range2 = PatchBouts.InPatchIndx(i-1,2):PatchBouts.InPatchIndx(i,1);
            if length(range2) > 300
                Behavior_C_leave{j}(c1,301:600) = Behavior(j,range2(1:300));
            else
                Behavior_C_leave{j}(c1,301:300+length(range2)) = Behavior(j,range2);
                Behavior_C_leave{j}(c1,600-(300-length(range2)):600) = NaN;
            end
            end
        end
    else
        c2 = c2 + 1;
        for j = 1:4
            range1 = PatchBouts.InPatchIndx(i-1,2):(PatchBouts.InPatchIndx(i,1)-1);
            if length(range1) > 300
                Behavior_C_return{j}(c2,1:300) = Behavior(j,range1(end-299:end));
            else
                Behavior_C_return{j}(c2,1:300-length(range1)) = NaN;
                Behavior_C_return{j}(c2,300-length(range1)+1:300) = Behavior(j,range1);
            end
            range2 = PatchBouts.InPatchIndx(i,1):PatchBouts.InPatchIndx(i,2);
            if length(range2) > 300
                Behavior_C_return{j}(c2,301:600) = Behavior(j,range2(1:300));
            else
                Behavior_C_return{j}(c2,301:300+length(range2)) = Behavior(j,range2);
                Behavior_C_return{j}(c2,600-(300-length(range2)):600) = NaN;
            end
        end
    end
end
% 
% Behavior_new_Patch_Avg = []; Behavior_new_Patch_STD = []; Behavior_new_Patch_Ns = [];
% Behavior_return_Patch_Avg = []; Behavior_return_Patch_STD = []; Behavior_return_Patch_Ns = [];
% Behavior_leave_Patch_Avg = []; Behavior_leave_Patch_STD = []; Behavior_leave_Patch_Ns = [];
% for j = 1:4
%     [Behavior_new_Patch_Avg(j,:), Behavior_new_Patch_STD(j,:), Behavior_new_Patch_Ns(j,:)] = av_fins(Behavior_C_new{j});
%     [Behavior_return_Patch_Avg(j,:), Behavior_return_Patch_STD(j,:), Behavior_return_Patch_Ns(j,:)] = av_fins(Behavior_C_return{j});
%     [Behavior_leave_Patch_Avg(j,:), Behavior_leave_Patch_STD(j,:), Behavior_leave_Patch_Ns(j,:)] = av_fins(Behavior_C_leave{j});
% end
% 
% x = [1:600];x = x/frameRate;
% 
% figure;hold on;
% plot(x,Behavior_new_Patch_Avg(1,:),'color',clr1); 
% plot(x,Behavior_new_Patch_Avg(2,:),'color',clr2);
% plot(x,Behavior_new_Patch_Avg(3,:),'color',clr3);
% 
% figure;hold on;
% plot(x,Behavior_return_Patch_Avg(1,:),'color',clr1); 
% plot(x,Behavior_return_Patch_Avg(2,:),'color',clr2);
% plot(x,Behavior_return_Patch_Avg(3,:),'color',clr3);
% 
% figure;hold on;
% plot(x,Behavior_leave_Patch_Avg(1,:),'color',clr1); 
% plot(x,Behavior_leave_Patch_Avg(2,:),'color',clr2);
% plot(x,Behavior_leave_Patch_Avg(3,:),'color',clr3);
% 
%%


%%
save([Worm ' PatchAnalysis_Games'],'Worm','nosetip','NumPatches','PatchResidency','NoseToPatchDistances','t_WL','Closeness',...
    'binaryImage','PatchCentroids','PatchImages','PatchPixels','PatchEdges','PatchMajAxis','PatchAreas',...
    'CoM_Vel','CoM_Vel_Roam','PatchBouts','InPatch_Residency','TotalFeedingTime','PercOfPatchEaten',...
    'PercOfPatchVisited','NumberOfFeedingsPerPatch','PercentOfReturnsFreshEncounters','PercentOfReturnsFreshEncounters_Cumulative','FeedingTimePerVisit','FeedingTimePerPatchVisit',...
    'PatchCoveragePerVisit','PatchCoveragePerPatchVisit','FeedingSearch_Feeding_Times','FeedingBreak_Feeding_Times','NewFeedingVisitTime',...
    'NewFeedingVisitCoverage','FirstFeedingVisitTime','FirstFeedingVisitCoverage','ReturnFeedingVisitTime','ReturnFeedingVisitCoverage',...
    'NonFeedingPatchEncounter_FoodEaten','NonFeedingPatchEncounter_CoM','FeedingPatchEncounter_CoM','AnteriorAngleAmplitude_C_new',...
    'AnteriorAngleAmplitude_C_return','AnteriorAngleAmplitude_C_leave','AnteriorAngleRatio_C_new','AnteriorAngleRatio_C_return',...
    'AnteriorAngleRatio_C_leave','Ant_Bend_Pos','Ant_Bend_Neg','PrimaryBehavior','BehaviorEdges','BehaviorEdges_Logical',...
    'BehaviorSwitches','BehSwitches_C_new','BehSwitches_C_return','BehSwitches_C_leave','DistFromPatch_C_new','DistFromPatch_C_return',...
    'DistFromPatch_C_leave','DistFromPatch_C_break','TimeToGetAway','TimeToReachDistanceFromPatch','TimeToGetAway_Break',...
    'LocalSearchToRoamCrossOver','LocalSearchToRoamSlopes','MaxDistanceForFeedingBreaks','MeanDistanceForFeedingBreaks','ApproachToPatch','ApproachToPatchFit',...
    'CoM_V_C_new','CoM_V_C_return','CoM_V_C_leave','CoM_V_C_break','CoM_Encounter_Feeding','CoM_Encounter_NonFeeding',...
    'CoM_Encounter_WithFood','CoM_Encounter_NoFood','SumAngles_C_new','SumAngles_C_return','SumAngles_C_leave','dAngles_C_new',...
    'dAngles_C_return','dAngles_C_leave','Behavior_C_new','Behavior_C_return','Behavior_C_leave','NoseToPatchCentroids','PercentOfPatchEatenWhenLeaving',...
    'DistFromCentroid_C_feeding','DistFromCentroid_C_leave','DistFromCentroid_C_leave_NumEvents','CoordsForLeavingEvents',...
    'ThetaForLeavingEvents','ThetaFromCentroid_leave','NoseToPatchCentroids_theta','LeavingEventsRadialRange','DistFromCentroidEvents');











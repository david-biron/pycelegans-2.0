% SegNum = 18; frameRate = 10; scrsz = get(0,'screensize'); 
% 
% TotSeg = 20;                    %how many segments to divide worm by
% SegNum = TotSeg - 2;            %number of segments in the analysis
% SegLngth = 100/TotSeg;          %the point length of each segment along worm spline.
% SEGMENTS = 1:SegNum;            %set the range of segments for the analysis
% FRAMES = 1:NumFrames;           %set the range to be analyzed
% 
% 
% %% Compute Sliding Window of Quiescence
% 
% QBoutDurSW = zeros(2,length(RT_L)) * NaN; %Orig,FullQuiesc,BodyQuiesc
% QBoutNumSW = zeros(2,length(RT_L));
% QBoutQFSW = zeros(2,length(RT_L))*NaN;
% 
% for i = 3001:length(RT_L)-3000
%     for j = 1:2
%         NaNfrac = sum(RT_L(i-3000:i+3000)==0)/length(RT_L(i-3000:i+3000));
%         % If 33% of the data is missing, don't compute
%         if NaNfrac < (1/3)
%             temp = QBoutDur(j,i-3000:i+3000);
%             QBoutDurSW(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
%             QBoutNumSW(j,i) = sum(QBoutNum(j,i-3000:i+3000));
%             temp = QboutAll_LRT(j,i-3000:i+3000);
%             tt = (RT_L(i-3000:i+3000)==1); temp = temp(tt); %take out points where no data
%             
%             QBoutQFSW(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
%         end
%     end
% end
% 
% QFSW = zeros(2,length(RT_L))*NaN;
% for i = 3001:length(RT_L)-3000
%     for j = 1:2
%         NaNfrac = sum(RT_L(i-3000:i+3000)==0)/length(RT_L(i-3000:i+3000));
%         if NaNfrac < (1/3)
%             temp = QuiescAll_LRT(j,i-3000:i+3000);
%             tt = (RT_L(i-3000:i+3000)==1); temp = temp(tt); %take out points where no data
%             QFSW(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
%         end
%     end
% end
% 
% Seg_QF = zeros(18,length(RT_L))*NaN;
% for i = 3001:length(RT_L)-3000
%     for k = 1:18
%         temp = QBraw_RT(k,i-3000:i+3000);
%         tt = (RT_L(i-3000:i+3000)==1); temp = temp(tt); %take out points where no data
%         Seg_QF(k,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
%     end
% end
% 
% NumSeg_QF = zeros(18,length(RT_L))*NaN;
% for i = 3001:length(RT_L)-3000
%     for k = 1:18
%         temp = QCollective(k,i-3000:i+3000);
%         tt = (RT_L(i-3000:i+3000)==1); temp = temp(tt); %take out points where no data
%         NumSeg_QF(k,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
%     end
% end
% 
% %% Sliding Window For Behavior
% 
% BehaviorPrc10min = zeros(4,length(RT_L))*NaN; %FRW,REV,DWELL,QUIESC
% FRWVelocities = zeros(3,length(RT_L))*NaN; %All,Full,Mix
% REVVelocities = zeros(3,length(RT_L))*NaN;
% FRWFrequencies = zeros(3,length(RT_L))*NaN;
% REVFrequencies = zeros(3,length(RT_L))*NaN;
% BendInitFrequency = zeros(1,length(RT_L))*NaN;
% BendDuration = zeros(1,length(RT_L))*NaN;
% 
% for i = 3001:length(RT_L)-3000
%     NaNfrac = sum(RT_L(i-3000:i+3000)==0)/length(RT_L(i-3000:i+3000));
%     % If 33% of the data is missing, don't compute
%     if NaNfrac < (1/3)
%         for j = 1:4
%             temp = Behavior(j,i-3000:i+3000);
%             BehaviorPrc10min(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
%         end
%         for j = 1:3
%             temp = FRWVelocityRT(j,i-3000:i+3000);
%             if sum(temp(isfinite(temp))) > 600
%                 FRWVelocities(j,i) = mean(temp(isfinite(temp)));
%             end
%             temp = REVVelocityRT(j,i-3000:i+3000);
%             if sum(temp(isfinite(temp))) > 300
%                 REVVelocities(j,i) = mean(temp(isfinite(temp)));
%             end
%             temp = FRWStarts(j,i-3000:i+3000);
%             FRWFrequencies(j,i) = sum(temp)/600;
%             temp = REVStarts(j,i-3000:i+3000);
%             REVFrequencies(j,i) = sum(temp)/600;
%         end
%         BendInitFrequency(i) = sum(BendStarts(i-3000:i+3000))/600;
%         temp = BendDur(i-3000:i+3000);
%         BendDuration(i) = mean(temp(isfinite(temp)));
%     end
% end
% 
% WaveCoherenceSW = zeros(2,length(RT_L))*NaN;
% for i = 3001:length(RT_L)-3000
%     NaNfrac = sum(RT_L(i-3000:i+3000)==0)/length(RT_L(i-3000:i+3000));
%     % If 33% of the data is missing, don't compute
%     if NaNfrac < (1/3)
%         for j = 1:2
%             temp = WaveCoherence(j,i-3000:i+3000);
%             WaveCoherenceSW(j,i) = mean(temp(isfinite(temp)));
%         end
%     end
% end
% 
% EdgePrcSW10min = NaN(1,length(RT_L));
% for i = 3001:length(RT_L)-3000
%     NaNfrac = sum(RT_L(i-3000:i+3000)==0)/length(RT_L(i-3000:i+3000));
%     % If 33% of the data is missing, don't compute
%     if NaNfrac < (1/3)
%         temp = EdgePrc(i-3000:i+3000);
%         EdgePrcSW10min(i) = mean(temp(isfinite(temp)));
%     end
% end
% 
% %% Worm Length
% WL_Slid = zeros(1,length(RT_L))*NaN;
% 
% for i = 3001:length(RT_L)
%     if isfinite(SlidWindow(1,i))
%         temp = WormLengths(SlidWindow(1,i):SlidWindow(2,i));
%         WL_Slid(i) = mean(temp(isfinite(temp)));
%     else
%         break;
%     end
% end

%% Sum of Angles -- general 'bendiness' of worm
% AngSumSlid = zeros(1,length(RT_L))*NaN;
AngSumAbsSlid = zeros(1,length(RT_L))*NaN;
for i = 3001:length(RT_L)
    if isfinite(SlidWindow(1,i))
%         AngSumSlid(i) = mean(AngSum(SlidWindow(1,i):SlidWindow(2,i)));
        AngSumAbsSlid(i) = mean(AngSumAbs(SlidWindow(1,i):SlidWindow(2,i)));
    else
        break;
    end
end

%% diff(angles) -- general activity of worm
%Calculate mean diff(angles) for each segment with a 10 min window
DAng = zeros(SegNum,length(RT))*NaN;
for j = SEGMENTS
    for i = 3001:length(RT_L)
        if isfinite(SlidWindow(1,i))
            DAng(j,i) = mean(dAngles(j,(SlidWindow(1,i):SlidWindow(2,i))));
        else
            break;
        end
    end
end

%Calculate the mean diff(ang) over all segments
MeanDAng = mean(DAng);

%% Center of Mass
CoM_AreaSW = zeros(1,length(RT_L))*NaN;
for i = 3001:length(RT_L)-3000
    x = CoM_Area(i-3000:i+3000);
    x = x(isfinite(x));

    if ~isempty(x) 
        CoM_AreaSW(i) = mean(x);
    end
end

%%
save([Worm ' BehaviorSlidWindow'],'SlidWindow','frameRate','scrsz','Worm','SegNum','FRAMES','SegLngth','NumFrames','FRAMES','SEGMENTS','frames','RT','RT_L','IN','OUT',...
    'QBoutDurSW','QBoutNumSW','QBoutQFSW','QFSW','Seg_QF','NumSeg_QF',...
    'BehaviorPrc10min','FRWVelocities','REVVelocities','FRWFrequencies','REVFrequencies','BendInitFrequency','BendDuration',...
    'WaveCoherenceSW','WL_Slid','AngSumAbsSlid',...
    'DAng','MeanDAng','CoM_AreaSW','EdgePrcSW10min');



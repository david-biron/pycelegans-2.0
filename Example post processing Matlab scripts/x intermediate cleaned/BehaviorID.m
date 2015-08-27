% getsize = @(x)x(1);
SegNum = 18; frameRate = 10; scrsz = get(0,'screensize'); 
TotSeg = 20;                    %how many segments to divide worm by
SegNum = TotSeg - 2;            %number of segments in the analysis
SegLngth = 100/TotSeg;          %the point length of each segment along worm spline.
SEGMENTS = 1:SegNum;            %set the range of segments for the analysis

%% Find regions of quiescence for each segment
Q = {};
Qedges = {};
for i = SEGMENTS
    Q{i} = (dAngles(i,:)<0.001);        %find where diff(angles) is <0.002
    Qedges{i} = (find(diff(Q{i})));     %determine the edges
end

% store the frame numbers of segment specific quiescence bouts
Qbouts = {};
for j = SEGMENTS
    k = 1;
    for i = 1:length(Qedges{j})-1
        if Q{j}(Qedges{j}(i)+1) == 1
            Qbouts{j}(k,:) = [(Qedges{j}(i)+1) (Qedges{j}(i+1)+1)];
            k = k + 1;
        end
    end
end

%Filter out very short Qbouts (just one frame)
QboutsF = {};
for j = SEGMENTS
    k = 1;
    for i = 1:length(Qbouts{j})
        if diff(Qbouts{j}(i,:)) > 1        %Length Filter for Qbouts
            QboutsF{j}(k,:) = Qbouts{j}(i,:);
            k = k + 1;
        end
    end
end

%Include  inactive disruptions in Qbouts (if max(dAng) < 0.3)
QboutsF2 = {};
for j = SEGMENTS
    k = 1;
    QboutsF2{j}(k,:) = QboutsF{j}(1,:);
    for i = 2:length(QboutsF{j})
        x = abs(max(AHA_S05{5}(j,QboutsF{j}(i-1,end):QboutsF{j}(i,1))) - min(AHA_S05{5}(j,QboutsF{j}(i-1,end):QboutsF{j}(i,1))));
        if x > 0.3
            k = k + 1;
            QboutsF2{j}(k,:) = QboutsF{j}(i,:);
        else
            QboutsF2{j}(k,2) = QboutsF{j}(i,2);
        end
    end
end

%Store a binary realtime structure with raw Quiescence by Segment data
QBraw_RT = zeros(18,length(RT_L));
for j = SEGMENTS
    for i = 1:size(QboutsF2{j},1)
        QBraw_RT(j,frames(QboutsF2{j}(i,1)+1):frames(QboutsF2{j}(i,2)+1)) = 1;
    end
end
QBraw_RT = logical(QBraw_RT);

CumSumSegQuiesc = [];
for j = 1:18
    CumSumSegQuiesc(j,:) = cumsum(QBraw_RT(j,:));
end

%Store how many segments in Quiesc
QCollective = zeros(18,length(RT_L));
for i = 1:length(RT_L)
    for k = 1:18
        if sum(QBraw_RT(1:end,i)) == k
            QCollective(k,i) = 1;
        else
            QCollective(k,i) = 0;
        end
    end
end

CumSumNumSegQuiesc = [];
for j = 1:18
    CumSumNumSegQuiesc(j,:) = cumsum(QCollective(j,:));
end
%%%%%%%%%%%%%%%
QboutsFilt = {};
for j = SEGMENTS
    k = 1;
    for i = 1:length(QboutsF2{j})
        if diff(QboutsF2{j}(i,:)) > 20        %Length Filter for Segment Qbouts (2 sec)
            QboutsFilt{j}(k,:) = QboutsF2{j}(i,:);
            k = k + 1;
        end
    end
end

%Construct a boolean array to filter Qbouts
Qbouts_L = zeros(SegNum,length(AHA_S05{1}));
Mbouts_L = ones(SegNum,length(AHA_S05{1}));
for j = SEGMENTS
    for i = 1:size(QboutsFilt{j},1)
        Qbouts_L(j,QboutsFilt{j}(i,1):QboutsFilt{j}(i,end)) = 1;
        Mbouts_L(j,QboutsFilt{j}(i,1):QboutsFilt{j}(i,end)) = 0;
    end
    Medges{j} = (find(diff(Mbouts_L(j,:)))); 
end
Qbouts_L = logical(Qbouts_L);
Mbouts_L = logical(Mbouts_L);

%Compute Mbouts
Mbouts = {};
for j = SEGMENTS
    k = 1;
    for i = 1:length(Medges{j})-1
        if Mbouts_L(j,Medges{j}(i)+1) == 1
            Mbouts{j}(k,:) = [(Medges{j}(i)+1) (Medges{j}(i+1)+1)];
            k = k + 1;
        end
    end
end

%real time vector indicating presence of Qbouts
QB_RT = zeros(j,length(RT_L));
for j = SEGMENTS
    for i = 1:size(QboutsFilt{j},1)
        QB_RT(j,frames(QboutsFilt{j}(i,1)+1):frames(QboutsFilt{j}(i,2)+1)) = 1;
    end
end
QB_RT = logical(QB_RT);
QB_RT_SumSeg = sum(QB_RT(:,:));

%% identify periods of worm Quiescence

%Store the position of allQbouts
%   1: Original definition (>16 segments including head; allowing 2 segments to move, excluding the head)
%   2: Body Quiescence > 12 segments (allowing 5 segments to move
QboutAll_L = zeros(2,length(AHA_S05{1}));
for i = 1:length(Qbouts_L(1,:))
    if sum(Qbouts_L(2:end,i)) > 15 && Qbouts_L(1,i) == 1
        QboutAll_L(1,i) = 1;
    else
        QboutAll_L(1,i) = 0;
    end
    if sum(Qbouts_L(1:end,i)) > 12
        QboutAll_L(2,i) = 1;
    else
        QboutAll_L(2,i) = 0;
    end
        
end
QboutAll_L = logical(QboutAll_L);
MboutAll_L = ~QboutAll_L;

%%%%%%%%%%%%%%%%
QuiescAll_LRT = zeros(2,length(RT_L));
for i = 1:length(RT_L)
    if sum(QBraw_RT(2:end,i)) > 15 && QBraw_RT(1,i) == 1
        QuiescAll_LRT(1,i) = 1;
    else
        QuiescAll_LRT(1,i) = 0;
    end
    if sum(QBraw_RT(1:end,i)) > 12
        QuiescAll_LRT(2,i) = 1;
    else
        QuiescAll_LRT(2,i) = 0;
    end    
end
QuiescAll_LRT = logical(QuiescAll_LRT);
MotAll_L = ~QuiescAll_LRT;
%%%%%%%%%%%%%%%%%%%%%%

% Pull out QboutAll instances
QboutAll_t = {};
Qall_edges = {};
for j = 1:2
    Qall_edges{j} = (find(diff(QboutAll_L(j,:))));
    k = 1;
    for i = 1:length(Qall_edges{j})-1
        if QboutAll_L(j,Qall_edges{j}(i)+1) == 1
            QboutAll_t{j}(k,:) = [(Qall_edges{j}(i)+1) (Qall_edges{j}(i+1)+1)];
            k = k + 1;
        end
    end
end

%Filter out short QboutsAll
QboutAll = {};
for j = 1:2
    k = 1;
    for i = 1:size(QboutAll_t{j}(:,:),1)
        if diff(QboutAll_t{j}(i,:)) > 20        %Length Filter for Qbouts
            QboutAll{j}(k,:) = QboutAll_t{j}(i,:);
            k = k + 1;
        end
    end
end

QboutAll_LRT = zeros(2,length(RT));
QboutAll_L = zeros(2,length(AHA_S05{1}));
for j = 1:2
    for i = 1:size(QboutAll{j}(:,:),1)
        QboutAll_LRT(j,frames(QboutAll{j}(i,1)+1):frames(QboutAll{j}(i,2)+1)) = 1;
        QboutAll_L(j,QboutAll{j}(i,1):QboutAll{j}(i,2)) = 1;
    end
end
QboutAll_LRT = logical(QboutAll_LRT);
QboutAll_L = logical(QboutAll_L);
MboutAll_LRT = ~QboutAll_LRT;
MboutAll_L = ~QboutAll_L;

MboutAll = {};
Mall_edges = {};
for j = 1:2
    Mall_edges{j} = (find(diff(MboutAll_L(j,:))));
    k = 1;
    for i = 1:length(Mall_edges{j})-1
        if MboutAll_L(j,Mall_edges{j}(i)+1) == 1
            MboutAll{j}(k,:) = [(Mall_edges{j}(i)+1) (Mall_edges{j}(i+1)+1)];
            k = k + 1;
        end
    end
end

clear('Q','Qedges','QT','QTedges','Qbouts','QboutsF','QboutsF2','Medges','Qall_edges','Mall_edges','QboutAll_t');

%% Store QBoutData in bout# referenced structure
QBoutData = {};
if ~isempty(QboutsFilt)
    for j = SEGMENTS
        for i = 1:size(QboutsFilt{j},1)
            QBoutData{j,i}.Ang = AHA_S05{1}(j,QboutsFilt{j}(i,1):QboutsFilt{j}(i,(end)));
            QBoutData{j,i}.absAng = abs(AHA_S05{1}(j,QboutsFilt{j}(i,1):QboutsFilt{j}(i,(end))));
            QBoutData{j,i}.T = frames((QboutsFilt{j}(i,1)+1):(QboutsFilt{j}(i,(end))+1));
            QBoutData{j,i}.Frame = QboutsFilt{j}(i,1):QboutsFilt{j}(i,(end));
            QBoutData{j,1}.Num = i;
        end
    end
end
QAllBoutData = {};
if ~isempty(QboutAll)
    for k = 1:2
        for j = SEGMENTS
            for i = 1:size(QboutAll{k},1)
                QAllBoutData{k,j,i}.Ang = AHA_S05{1}(j,QboutAll{k}(i,1):QboutAll{k}(i,2));
                QAllBoutData{k,j,i}.absAng = abs(AHA_S05{1}(j,QboutAll{k}(i,1):QboutAll{k}(i,2)));
                QAllBoutData{k,j,i}.T = frames(QboutAll{k}(i,1)+1:QboutAll{k}(i,2)+1);
                QAllBoutData{k,j,i}.Frame = QboutAll{k}(i,1):QboutAll{k}(i,2);
                QAllBoutData{k,j,1}.Num = i;
            end
        end
    end
end

MBoutData = {};
if ~isempty(Mbouts)
    for j = SEGMENTS
        for i = 1:size(Mbouts{j},1)
            MBoutData{j,i}.Ang = AHA_S05{1}(j,Mbouts{j}(i,1):Mbouts{j}(i,(end)));
            MBoutData{j,i}.absAng = abs(AHA_S05{1}(j,Mbouts{j}(i,1):Mbouts{j}(i,(end))));
            MBoutData{j,i}.T = frames(Mbouts{j}(i,1)+1:Mbouts{j}(i,(end))+1);
            MBoutData{j,i}.Frame = Mbouts{j}(i,1):Mbouts{j}(i,(end));
            MBoutData{j,1}.Num = i;
        end
    end
end
MAllBoutData = {};
if ~isempty(MboutAll)
    for k = 1:2
        for j = SEGMENTS
            for i = 1:size(MboutAll{k},1)
                MAllBoutData{k,j,i}.Ang = AHA_S05{1}(j,MboutAll{k}(i,1):MboutAll{k}(i,2));
                MAllBoutData{k,j,i}.absAng = abs(AHA_S05{1}(j,MboutAll{k}(i,1):MboutAll{k}(i,2)));
                MAllBoutData{k,j,i}.T = frames(MboutAll{k}(i,1)+1:MboutAll{k}(i,2)+1);
                MAllBoutData{k,j,i}.Frame = MboutAll{k}(i,1):MboutAll{k}(i,2);
                MAllBoutData{k,j,1}.Num = i;
            end
        end
    end
end

%% Compile QBnextMB and MBnextQB
QB = {};
MB = {};
if ~isempty(QAllBoutData{1})
    s = size(QAllBoutData); s = s(3);
    for k = 1:2
        j = 1;
        for i = 1:s
            if ~isempty(QAllBoutData{k,1,i})
                QB{k}(j,:) = [QAllBoutData{k,1,i}.T(1) ((QAllBoutData{k,1,i}.T(end)-QAllBoutData{k,1,i}.T(1))/frameRate)];
                j = j + 1;
            end
        end
    end
    s = size(MAllBoutData); s = s(3);
    for k = 1:2
        j = 1;
        for i = 1:s
            if ~isempty(MAllBoutData{k,1,i})
                MB{k}(j,:) = [(MAllBoutData{k,1,i}.T(1)) ((MAllBoutData{k,1,i}.T(end)-MAllBoutData{k,1,i}.T(1))/frameRate)];
                j = j + 1;
            end
        end
    end
end
QBnextMB = {}; %for each QuiescID: {1-3}, [duration of QB, duration of next MB, start time of QB, start time of MB]
MBnextQB = {};

if ~isempty(QB)
    for k = 1:2
        if ~isempty(QB{k}) && ~isempty(MB{k}) && ~isempty(QB{k})
            if QB{k}(1,1) < MB{k}(1,1)
                for i = 1:min([size(QB{k},1) size(MB{k},1)])
                    QBnextMB{k}(i,:) = [QB{k}(i,2) MB{k}(i,2) QB{k}(i,1) MB{k}(i,1)];
                end
                for i = 1:min([size(QB{k},1) size(MB{k},1)])-1
                    MBnextQB{k}(i,:) = [MB{k}(i,2) QB{k}(i+1,2) MB{k}(i,1) QB{k}(i,1)];
                end
            end
            if QB{k}(1,1) > MB{k}(1,1)
                for i = 1:min([size(QB{k},1) size(MB{k},1)])
                    MBnextQB{k}(i,:) = [MB{k}(i,2) QB{k}(i,2) MB{k}(i,1) QB{k}(i,1)];
                end
                for i = 1:min([size(QB{k},1) size(MB{k},1)])-1
                    QBnextMB{k}(i,:) = [QB{k}(i,2) MB{k}(i+1,2) QB{k}(i,1) MB{k}(i,1)];
                end
            end
        end
    end
end
clear('QB','MB','s');
%%
QBoutDur = zeros(2,length(RT_L)) * NaN; %Original, All, Body
QBoutNum = zeros(2,length(RT_L));

if ~isempty(QAllBoutData{1})
    s = size(QAllBoutData); s = s(3);
    for k = 1:2
        for i = 1:s
            if ~isempty(QAllBoutData{k,1,i})
                QBoutDur(k,QAllBoutData{k,1,i}.T(1)) = ((QAllBoutData{k,1,i}.T(end)-QAllBoutData{k,1,i}.T(1))/frameRate);
                QBoutNum(k,QAllBoutData{k,1,i}.T(1)) = 1;
            end
        end
    end
end

%% Identify Peaks in body posture

PeaksW = {};
PeaksWIndx = {};
PeaksWNeg = {};
PeaksWNegIndx = {};
for i = 1:length(RT_L)
    
    if isempty(find(frames==i-1))
        PeaksW{i} = [];
        PeaksWIndx{i} = [];
        PeaksWNeg{i} = [];
        PeaksWNegIndx{i} = [];
    else
        x1 = find(frames==i-1);
        yy = [];
        xx = [];
        for j = SEGMENTS
            yy = [yy AHA_S05{3}(j,x1)];
            xx = [xx j];
        end
        yyNeg = yy*-1;
        
        yy = smooth(yy,3);
        yyNeg = smooth(yyNeg,3);
        
        [P Pindx] = findpeaks(yy,'minpeakdistance',3);
        PeaksW{i} = P';
        PeaksWIndx{i} = (Pindx);
        
        [PNeg PNegindx] = findpeaks(yyNeg,'minpeakdistance',3);
        PeaksWNeg{i} = PNeg';
        PeaksWNegIndx{i} = (PNegindx);
    end
end

%Pull out positive Waves
k = 0;
PWaves = {};
cc = 1;
Wid = 0;
WidT = zeros(1,18);
for i = 1:length(RT_L)-11
    if i >= max(frames)
        break;
    end
    if isempty(PeaksWIndx{i})
        k = k + 1;
        continue;
    end
    if k > 10
        WidT = zeros(1,18);
        k = 0;
    else
        k = 0;
    end
    
    tt = 0;
    if isempty(PeaksWIndx{i+1})
        for ii = 2:10
            if ~isempty(PeaksWIndx{i+ii})
                tt = 1;
                break;
            end
        end
    else
        tt = 1;
        ii = 1;
    end
    if tt == 0
        continue;
    end
    
    WJc = 1;
    WJ = [];
    for j = 1:length(PeaksWIndx{i})
        PeakSeg = PeaksWIndx{i}(j);
        PeakAng = PeaksW{i}(j);

        %check to see if there is an already started wave
        if WidT(PeakSeg) == 0
            %establish new Wave# pointer
            Wid = cc;
            cc = cc + 1;
            
            WidT(PeakSeg) = Wid;
            
            PWaves{(Wid)}.A = [PeakAng];
            PWaves{(Wid)}.T = [i];
            PWaves{(Wid)}.S = [PeakSeg];
        else
            %If theres already a wave there, set up pointer to its ID
            Wid = WidT(PeakSeg);
        end
        
        %Determine brackets to help look for next Segment for the wave
        
        %determine brackets based on presence of opposite waves in vicinity
        x2 = PeaksWNegIndx{i}-PeakSeg;
        if isempty(x2)
            bracket1a = 1;
            bracket2a = 18;
        else
            if isempty(max(x2(x2<0)))
                bracket1a = 1;
                bracket2a = PeaksWNegIndx{i}(1);
            else
                if isempty(min(x2(x2>0)))
                    bracket1a = PeaksWNegIndx{i}(end);
                    bracket2a = 18;
                else
                    bracket1a = PeaksWNegIndx{i}(find(x2 == max(x2(x2<0))));
                    bracket2a = PeaksWNegIndx{i}(find(x2 == min(x2(x2>0))));
                end
            end
        end
        %Do the same for the next frame to be careful
        x2 = PeaksWNegIndx{i+ii}-PeakSeg;
        if isempty(x2)
            bracket1b = 1;
            bracket2b = 18;
        else
            if isempty(max(x2(x2<0)))
                bracket1b = 1;
                bracket2b = PeaksWNegIndx{i+ii}(1);
            else
                if isempty(min(x2(x2>0)))
                    bracket1b = PeaksWNegIndx{i+ii}(end);
                    bracket2b = 18;
                else
                    bracket1b = PeaksWNegIndx{i+ii}(find(x2 == max(x2(x2<0))));
                    bracket2b = PeaksWNegIndx{i+ii}(find(x2 == min(x2(x2>0))));
                end
            end
        end
        
        bracket1 = max(bracket1a,bracket1b);
        bracket2 = min(bracket2a,bracket2b);
        
        %look over all peaks in the next frame
        ff = 0;
        fff = [];
        for jj = 1:length(PeaksWIndx{i+ii})
            if (PeaksWIndx{i+ii}(jj) < bracket2) && (PeaksWIndx{i+ii}(jj) > bracket1)
                ff = ff + 1;
                fff(ff) = jj;
            end
        end
        
        if ff > 1   %if forks, pick the closest offpath.
            dff = abs(fff-PeakSeg);
            jjj = fff(find(dff==min(dff)));
            WJ(WJc,1) = Wid;
            WJ(WJc,2) = jjj;
            WJc = WJc + 1;
        else
            if ff == 1
                jjj = fff(ff);
                WJ(WJc,1) = Wid;
                WJ(WJc,2) = jjj;
                WJc = WJc + 1;
            else
                WidT(find(WidT==Wid)) = 0;
            end
        end
    end
    
    if ~isempty(WJ)
        for kk = 1:length(WJ(:,1));
            
            WidT(find(WidT==WJ(kk,1))) = 0;
            WidT(PeaksWIndx{i+ii}(WJ(kk,2))) = WJ(kk,1);
            
            PWaves{(WJ(kk,1))}.A = [PWaves{(WJ(kk,1))}.A PeaksW{i+ii}(WJ(kk,2))];
            PWaves{(WJ(kk,1))}.T = [PWaves{(WJ(kk,1))}.T (i+ii)];
            PWaves{(WJ(kk,1))}.S = [PWaves{(WJ(kk,1))}.S PeaksWIndx{i+ii}(WJ(kk,2))];
            
        end
    end
end
k = 1;
Waves1 = {};
for i = 1:length(PWaves)
    if length(PWaves{i}.S) > 3 && abs(max(PWaves{i}.S)-min(PWaves{i}.S)) > 3
        Waves1{k} = PWaves{i};
        k = k + 1;
    end
end

%Pull out negative waves
k = 0;
NWaves = {};
cc = 1;
Wid = 0;
WidT = zeros(1,18);
for i = 1:length(RT_L)-11
    if i >= max(frames)
        break;
    end
    if isempty(PeaksWNegIndx{i})
        k = k + 1;
        continue;
    end
    if k > 10
        WidT = zeros(1,18);
        k = 0;
    else
        k = 0;
    end
    
    tt = 0;
    if isempty(PeaksWNegIndx{i+1})
        for ii = 2:10
            if ~isempty(PeaksWNegIndx{i+ii})
                tt = 1;
                break;
            end
        end
    else
        tt = 1;
        ii = 1;
    end
    if tt == 0
        continue;
    end
    
    WJc = 1;
    WJ = [];
    for j = 1:length(PeaksWNegIndx{i})
        PeakSeg = PeaksWNegIndx{i}(j);
        PeakAng = PeaksWNeg{i}(j);

        %check to see if there is an already started wave
        if WidT(PeakSeg) == 0
            %establish new Wave# pointer
            Wid = cc;
            cc = cc + 1;
            
            WidT(PeakSeg) = Wid;
            
            NWaves{(Wid)}.A = [PeakAng];
            NWaves{(Wid)}.T = [i];
            NWaves{(Wid)}.S = [PeakSeg];
        else
            %If theres already a wave ther, set up pointer to its ID
            Wid = WidT(PeakSeg);
        end
        
        %Determine brackets to help look for next Segment for the wave
        
        %determine brackets based on presence of opposite waves in vicinity
        x2 = PeaksWIndx{i}-PeakSeg;
        if isempty(x2)
            bracket1a = 1;
            bracket2a = 18;
        else
            if isempty(max(x2(x2<0)))
                bracket1a = 1;
                bracket2a = PeaksWIndx{i}(1);
            else
                if isempty(min(x2(x2>0)))
                    bracket1a = PeaksWIndx{i}(end);
                    bracket2a = 18;
                else
                    bracket1a = PeaksWIndx{i}(find(x2 == max(x2(x2<0))));
                    bracket2a = PeaksWIndx{i}(find(x2 == min(x2(x2>0))));
                end
            end
        end
        %Do the same for the next frame to be careful
        x2 = PeaksWIndx{i+ii}-PeakSeg;
        if isempty(x2)
            bracket1b = 1;
            bracket2b = 18;
        else
            if isempty(max(x2(x2<0)))
                bracket1b = 1;
                bracket2b = PeaksWIndx{i+ii}(1);
            else
                if isempty(min(x2(x2>0)))
                    bracket1b = PeaksWIndx{i+ii}(end);
                    bracket2b = 18;
                else
                    bracket1b = PeaksWIndx{i+ii}(find(x2 == max(x2(x2<0))));
                    bracket2b = PeaksWIndx{i+ii}(find(x2 == min(x2(x2>0))));
                end
            end
        end
        
        bracket1 = max(bracket1a,bracket1b);
        bracket2 = min(bracket2a,bracket2b);
        
        %look over all peaks in the next frame
        ff = 0;
        fff = [];
        for jj = 1:length(PeaksWNegIndx{i+ii})
            if (PeaksWNegIndx{i+ii}(jj) < bracket2) && (PeaksWNegIndx{i+ii}(jj) > bracket1)
                ff = ff + 1;
                fff(ff) = jj;
            end
        end
        
        if ff > 1   %if forks, pick the closest offpath.
            dff = abs(fff-PeakSeg);
            jjj = fff(find(dff==min(dff)));
            WJ(WJc,1) = Wid;
            WJ(WJc,2) = jjj;
            WJc = WJc + 1;
        else
            if ff == 1
                jjj = fff(ff);
                WJ(WJc,1) = Wid;
                WJ(WJc,2) = jjj;
                WJc = WJc + 1;
            else
                WidT(find(WidT==Wid)) = 0;
            end
        end
    end
    
    if ~isempty(WJ)
        for kk = 1:length(WJ(:,1));
            
            WidT(find(WidT==WJ(kk,1))) = 0;
            WidT(PeaksWNegIndx{i+ii}(WJ(kk,2))) = WJ(kk,1);
            
            NWaves{(WJ(kk,1))}.A = [NWaves{(WJ(kk,1))}.A PeaksWNeg{i+ii}(WJ(kk,2))];
            NWaves{(WJ(kk,1))}.T = [NWaves{(WJ(kk,1))}.T (i+ii)];
            NWaves{(WJ(kk,1))}.S = [NWaves{(WJ(kk,1))}.S PeaksWNegIndx{i+ii}(WJ(kk,2))];
            
        end
    end
end
k = 1;
Waves2 = {};
for i = 1:length(NWaves)
    if length(NWaves{i}.S) > 3 && abs(max(NWaves{i}.S)-min(NWaves{i}.S)) > 3
        Waves2{k} = NWaves{i};
        k = k + 1;
    end
end

% Clean Up waves and pull out directionality
for i = 1:length(Waves1)
    xx = {};
    W = Waves1{i}.S;
    T = Waves1{i}.T;
    I = (1:length(Waves1{i}.T));
    Q = QboutAll_LRT(2,T);
    Quiesc = QboutAll_LRT(1,T);
    Waves1{i}.QT = T(Q);
    Waves1{i}.QS = W(Q);
    Waves1{i}.QI = I(Q);
    Waves1{i}.QuiescT = T(Quiesc);
    Waves1{i}.QuiescS = W(Quiesc);
    Waves1{i}.QuiescI = I(Quiesc);
    
    Qdiff = diff(Q);
    
    Wdiff = diff(W);
    Wdiff(Wdiff > 0) = 1;
    Wdiff(Wdiff < 0) = -1;
    Wdiff(Qdiff~=0) = 20;
%     W = [W(Wdiff~=0) W(end)];
%     T = [T(Wdiff~=0) T(end)];
%     I = [I(Wdiff~=0) I(end)];
    I1 = [I(Wdiff~=0)];
    I2 = [I(find(Wdiff~=0)+1)];
    
    Wdiff = (Wdiff(Wdiff~=0));
    
    Edge_Wdiff = find(diff(Wdiff));
    
    if isempty(Edge_Wdiff)
        xx{1}.S = W([I1 I2(end)]);
        xx{1}.T = T([I1 I2(end)]);
        xx{1}.I = I([I1 I2(end)]);
    else
        xx{1}.S = W([I1(1:Edge_Wdiff(1)) I2(Edge_Wdiff(1))]);
        xx{1}.T = T([I1(1:Edge_Wdiff(1)) I2(Edge_Wdiff(1))]);
        xx{1}.I = I([I1(1:Edge_Wdiff(1)) I2(Edge_Wdiff(1))]);
        if length(Edge_Wdiff) >= 2
            for j = 2:length(Edge_Wdiff)
                xx{j}.S = W([I1(1+Edge_Wdiff(j-1):Edge_Wdiff(j)) I2(Edge_Wdiff(j))]);
                xx{j}.T = T([I1(1+Edge_Wdiff(j-1):Edge_Wdiff(j)) I2(Edge_Wdiff(j))]);
                xx{j}.I = I([I1(1+Edge_Wdiff(j-1):Edge_Wdiff(j)) I2(Edge_Wdiff(j))]);
            end
        else
            j = 1;
        end
        xx{j+1}.S = W([I1(1+Edge_Wdiff(j):end) I2((end))]);
        xx{j+1}.T = T([I1(1+Edge_Wdiff(j):end) I2((end))]);
        xx{j+1}.I = I([I1(1+Edge_Wdiff(j):end) I2((end))]);
    end
    filt_xx = {};
    k = 1;
    for j = 1:length(xx)
        if length(xx{j}.S) > 3
            filt_xx{k} = xx{j};
            k = k + 1;
        end
    end
    PosteriorT = [];
    AnteriorT = [];
    PosteriorS = [];
    AnteriorS = [];
    ForwardEvents = [];
    k = 1;
    ReverseEvents = [];
    z = 1;
    for j = 1:length(filt_xx)
        if mean(diff(filt_xx{j}.S)) > 0 && (sum(QboutAll_LRT(2,filt_xx{j}.T)) == 0)
            PosteriorT = [PosteriorT filt_xx{j}.T];
            PosteriorS = [PosteriorS filt_xx{j}.S];
            ForwardEvents(k,:) = [filt_xx{j}.I(1) filt_xx{j}.I(end)];
            k = k + 1;
        else
            if mean(diff(filt_xx{j}.S)) < 0 && (sum(QboutAll_LRT(2,filt_xx{j}.T)) == 0)
                AnteriorT = [AnteriorT filt_xx{j}.T];
                AnteriorS = [AnteriorS filt_xx{j}.S];
                ReverseEvents(z,:) = [filt_xx{j}.I(1) filt_xx{j}.I(end)];
                z = z + 1;
            end
        end
    end
    Waves1{i}.FRWI = ForwardEvents;
    Waves1{i}.REVI = ReverseEvents;
    Waves1{i}.FRWT = PosteriorT;
    Waves1{i}.FRWS = PosteriorS;
    Waves1{i}.REVT = AnteriorT;
    Waves1{i}.REVS = AnteriorS;
    Waves1{i}.FRWF = length(PosteriorT)/(length(PosteriorT) + length(AnteriorT));
    Waves1{i}.REVF = length(AnteriorT)/(length(PosteriorT) + length(AnteriorT));
end

for i = 1:length(Waves2)
    xx = {};
    W = Waves2{i}.S;
    T = Waves2{i}.T;
    I = (1:length(Waves2{i}.T));
    Q = QboutAll_LRT(2,T);
    Quiesc = QboutAll_LRT(1,T);
    Waves2{i}.QT = T(Q);
    Waves2{i}.QS = W(Q);
    Waves2{i}.QI = I(Q);
    Waves2{i}.QuiescT = T(Quiesc);
    Waves2{i}.QuiescS = W(Quiesc);
    Waves2{i}.QuiescI = I(Quiesc);
  
    Qdiff = diff(Q);
    
    Wdiff = diff(W);
    Wdiff(Wdiff > 0) = 1;
    Wdiff(Wdiff < 0) = -1;
    Wdiff(Qdiff~=0) = 20;
    
    I1 = [I(Wdiff~=0)];
    I2 = [I(find(Wdiff~=0)+1)];
    
    Wdiff = (Wdiff(Wdiff~=0));
    
    Edge_Wdiff = find(diff(Wdiff));
    
    if isempty(Edge_Wdiff)
        xx{1}.S = W([I1 I2(end)]);
        xx{1}.T = T([I1 I2(end)]);
        xx{1}.I = I([I1 I2(end)]);
    else
        xx{1}.S = W([I1(1:Edge_Wdiff(1)) I2(Edge_Wdiff(1))]);
        xx{1}.T = T([I1(1:Edge_Wdiff(1)) I2(Edge_Wdiff(1))]);
        xx{1}.I = I([I1(1:Edge_Wdiff(1)) I2(Edge_Wdiff(1))]);
        if length(Edge_Wdiff) >= 2
            for j = 2:length(Edge_Wdiff)
                xx{j}.S = W([I1(1+Edge_Wdiff(j-1):Edge_Wdiff(j)) I2(Edge_Wdiff(j))]);
                xx{j}.T = T([I1(1+Edge_Wdiff(j-1):Edge_Wdiff(j)) I2(Edge_Wdiff(j))]);
                xx{j}.I = I([I1(1+Edge_Wdiff(j-1):Edge_Wdiff(j)) I2(Edge_Wdiff(j))]);
            end
        else
            j = 1;
        end
        xx{j+1}.S = W([I1(1+Edge_Wdiff(j):end) I2((end))]);
        xx{j+1}.T = T([I1(1+Edge_Wdiff(j):end) I2((end))]);
        xx{j+1}.I = I([I1(1+Edge_Wdiff(j):end) I2((end))]);
    end
    filt_xx = {};
    k = 1;
    for j = 1:length(xx)
        if length(xx{j}.S) > 3
            filt_xx{k} = xx{j};
            k = k + 1;
        end
    end
    PosteriorT = [];
    AnteriorT = [];
    PosteriorS = [];
    AnteriorS = [];
    ForwardEvents = [];
    k = 1;
    ReverseEvents = [];
    z = 1;
    for j = 1:length(filt_xx)
        if mean(diff(filt_xx{j}.S)) > 0 && (sum(QboutAll_LRT(2,filt_xx{j}.T)) == 0)
            PosteriorT = [PosteriorT filt_xx{j}.T];
            PosteriorS = [PosteriorS filt_xx{j}.S];
            ForwardEvents(k,:) = [filt_xx{j}.I(1) filt_xx{j}.I(end)];
            k = k + 1;
        else
            if mean(diff(filt_xx{j}.S)) < 0 && (sum(QboutAll_LRT(2,filt_xx{j}.T)) == 0)
                AnteriorT = [AnteriorT filt_xx{j}.T];
                AnteriorS = [AnteriorS filt_xx{j}.S];
                ReverseEvents(z,:) = [filt_xx{j}.I(1) filt_xx{j}.I(end)];
                z = z + 1;
            end
        end
    end
    Waves2{i}.FRWI = ForwardEvents;
    Waves2{i}.REVI = ReverseEvents;
    Waves2{i}.FRWT = PosteriorT;
    Waves2{i}.FRWS = PosteriorS;
    Waves2{i}.REVT = AnteriorT;
    Waves2{i}.REVS = AnteriorS;
    Waves2{i}.FRWF = length(PosteriorT)/(length(PosteriorT) + length(AnteriorT));
    Waves2{i}.REVF = length(AnteriorT)/(length(PosteriorT) + length(AnteriorT));
end

clear('P','Pindx','PNeg','PeakAng','PeakSeg','PNegindx','PWaves','NWaves','PeaksW','PeaksWIndx','PeaksWNeg','PeaksWNegIndx','Wid','WidT','AnteriorS','AnteriorT','PosteriorS','PosteriorT','Edge_Wdiff','ForwardEvents','ReverseEvents','Wdiff','Qdiff','WJ','WJc');
clear('bracket1','bracket1a','bracket1b','bracket2','bracket2a','bracket2b','cc')
clear('dff','ff','fff','jjj','filt_xx','xx','yy','yyNeg')
clear('Quiesc','W','T','I','Q');
%% Calculations involving Waves
for i = 1:length(Waves1)
    if ~isempty(Waves1{i}.FRWI)
        for j = 1:length(Waves1{i}.FRWI(:,1))
            Range = [Waves1{i}.FRWI(j,1) Waves1{i}.FRWI(j,2)];
            Waves1{i}.FRWV(j) = (sum(diff(Waves1{i}.S(Range(1):Range(2))))*5) / ((Waves1{i}.T(Range(2))-Waves1{i}.T(Range(1)))/10);
            Waves1{i}.FrwSegCov(j) = max(Waves1{i}.S(Range(1):Range(2))) - min(Waves1{i}.S(Range(1):Range(2)));
        end
    else
        Waves1{i}.FRWV = NaN;
        Waves1{i}.FrwSegCov = NaN;
    end
    if ~isempty(Waves1{i}.REVI)
        for j = 1:length(Waves1{i}.REVI(:,1))
            Range = [Waves1{i}.REVI(j,1) Waves1{i}.REVI(j,2)];
            Waves1{i}.REVV(j) = -1*(sum(diff(Waves1{i}.S(Range(1):Range(2))))*5) / ((Waves1{i}.T(Range(2))-Waves1{i}.T(Range(1)))/10);
            Waves1{i}.RevSegCov(j) = max(Waves1{i}.S(Range(1):Range(2))) - min(Waves1{i}.S(Range(1):Range(2)));
        end
    else
        Waves1{i}.REVV = NaN;
        Waves1{i}.RevSegCov = NaN;
    end
    Waves1{i}.Dur = (Waves1{i}.T(end)-Waves1{i}.T(1))/frameRate;
    Waves1{i}.EdgePrc = sum(OUT(Waves1{i}.T))/(sum(IN(Waves1{i}.T))+sum(OUT(Waves1{i}.T)));
end
for i = 1:length(Waves2)
    if ~isempty(Waves2{i}.FRWI)
        for j = 1:length(Waves2{i}.FRWI(:,1))
            Range = [Waves2{i}.FRWI(j,1) Waves2{i}.FRWI(j,2)];
            Waves2{i}.FRWV(j) = (sum(diff(Waves2{i}.S(Range(1):Range(2))))*5) / ((Waves2{i}.T(Range(2))-Waves2{i}.T(Range(1)))/10);
            Waves2{i}.FrwSegCov(j) = max(Waves2{i}.S(Range(1):Range(2))) - min(Waves2{i}.S(Range(1):Range(2)));
        end
    else
        Waves2{i}.FRWV = NaN;
        Waves2{i}.FrwSegCov = NaN;
    end
    if ~isempty(Waves2{i}.REVI)
        for j = 1:length(Waves2{i}.REVI(:,1))
            Range = [Waves2{i}.REVI(j,1) Waves2{i}.REVI(j,2)];
            Waves2{i}.REVV(j) = -1*(sum(diff(Waves2{i}.S(Range(1):Range(2))))*5) / ((Waves2{i}.T(Range(2))-Waves2{i}.T(Range(1)))/10);
            Waves2{i}.RevSegCov(j) = max(Waves2{i}.S(Range(1):Range(2))) - min(Waves2{i}.S(Range(1):Range(2)));
        end
    else
        Waves2{i}.REVV = NaN;
        Waves2{i}.RevSegCov = NaN;
    end
    Waves2{i}.Dur = (Waves2{i}.T(end)-Waves2{i}.T(1))/frameRate;
    Waves2{i}.EdgePrc = sum(OUT(Waves2{i}.T))/(sum(IN(Waves2{i}.T))+sum(OUT(Waves2{i}.T)));
end
%% Construct RT Pointer to WaveID

Waves_Pointer = cell(length(RT),2);

for i = 1:length(Waves1)
    for j = 1:length(Waves1{i}.T)
        Waves_Pointer{Waves1{i}.T(j),1} = [Waves_Pointer{Waves1{i}.T(j),1}; [i j]];
    end
end
for i = 1:length(Waves2)
    for j = 1:length(Waves2{i}.T)
        Waves_Pointer{Waves2{i}.T(j),2} = [Waves_Pointer{Waves2{i}.T(j),2}; [i j]];
    end
end

%% Pull out primary locomotive behavior in RT 
F = zeros(1,length(RT)); %Number of waves moving Forward along RT
R = zeros(1,length(RT));
D = zeros(1,length(RT));
Q = zeros(1,length(RT));
Unc = zeros(1,length(RT));
Behavior = zeros(4,length(RT)); %FRW,REV,DWELL,QUIESC
EdgePrc = NaN(1,length(RT));

for i = 1:length(RT)
    temp1 = Waves_Pointer{i,1};
    temp2 = Waves_Pointer{i,2};
    noWaveFlag1 = 0;
    noWaveFlag2 = 0;
    EdgePrcTemp = [];
    if size(temp1,1) ~= 0
        for j = 1:size(temp1,1)
            flag = 0;
            if ismember(temp1(j,2),Waves1{temp1(j,1)}.QuiescI)
%                 Q(i) = Q(i) + 1;
            else
                if ~isempty(Waves1{temp1(j,1)}.FRWI)
                    for k = 1:size(Waves1{temp1(j,1)}.FRWI,1)
                        if temp1(j,2) >= Waves1{temp1(j,1)}.FRWI(k,1) && temp1(j,2) <= Waves1{temp1(j,1)}.FRWI(k,2)
                            F(i) = F(i) + 1;
                            flag = 1;
                        end
                    end
                end
                if ~isempty(Waves1{temp1(j,1)}.REVI)
                    for k = 1:size(Waves1{temp1(j,1)}.REVI,1)
                        if temp1(j,2) >= Waves1{temp1(j,1)}.REVI(k,1) && temp1(j,2) <= Waves1{temp1(j,1)}.REVI(k,2)
                            R(i) = R(i) + 1;
                            flag = 1;
                        end
                    end
                end
                if flag == 0
                    D(i) = D(i) + 1;
                end
            end
            EdgePrcTemp = [EdgePrcTemp Waves1{temp1(j,1)}.EdgePrc];
        end
    else
        noWaveFlag1 = 1;
    end
    
    if size(temp2,1) ~= 0
        for j = 1:size(temp2,1)
            flag = 0;
            if ismember(temp2(j,2),Waves2{temp2(j,1)}.QuiescI)
%                 Q(i) = Q(i) + 1;
            else
                if ~isempty(Waves2{temp2(j,1)}.FRWI)
                    for k = 1:size(Waves2{temp2(j,1)}.FRWI,1)
                        if temp2(j,2) >= Waves2{temp2(j,1)}.FRWI(k,1) && temp2(j,2) <= Waves2{temp2(j,1)}.FRWI(k,2)
                            F(i) = F(i) + 1;
                            flag = 1;
                        end
                    end
                end
                if ~isempty(Waves2{temp2(j,1)}.REVI)
                    for k = 1:size(Waves2{temp2(j,1)}.REVI,1)
                        if temp2(j,2) >= Waves2{temp2(j,1)}.REVI(k,1) && temp2(j,2) <= Waves2{temp2(j,1)}.REVI(k,2)
                            R(i) = R(i) + 1;
                            flag = 1;
                        end
                    end
                end
                if flag == 0;
                    D(i) = D(i) + 1;
                end
            end
            EdgePrcTemp = [EdgePrcTemp Waves2{temp2(j,1)}.EdgePrc];
        end
    else
        noWaveFlag2 = 1;
    end
    
    if QboutAll_LRT(1,i) == 1
        Q(i) = size(temp1,1) + size(temp2,1);
        if Q(i) == 0
            Q(i) = 20;
        end
    end
    
    if F(i) ~= 0 && R(i) ~= 0
        Unc(i) = 1;
    end
    
    
    if RT_L(i) == 0
        F(i) = NaN;
        R(i) = NaN;
        D(i) = NaN;
        Q(i) = NaN;
        Unc(i) = NaN;
        Behavior(:,i) = NaN;
    end
    
    b = [F(i) R(i) D(i) Q(i)];
    
    bb = [];
    if b(4) ~= 0 && isfinite(b(4))
        bb = 4;
    else
        if b(1) ~= 0 && b(1)>b(2)%b(2) == 0
            bb = 1;
        else
            if b(2) ~= 0 && b(2)>b(1)%b(1) == 0
                bb = 2;
            else
                if isfinite(b(3))
                    bb = 3;
                end
            end
        end
    end

    if ~isempty(bb)
        Behavior(bb,i) = 1;
    end
    EdgePrc(i) = mean(EdgePrcTemp);
end
clear('bb','noWaveFlag1','noWaveFlag2','b','temp1','temp2','flag')
%% Construct RT vector with FRW and REV parameters and general bend properties

FRWVelocityRT = zeros(3,length(RT))*NaN; %All, Full, Mix
REVVelocityRT = zeros(3,length(RT))*NaN;
FRWEdgeFiltVel = zeros(6,length(RT))*NaN; %All_IN, All_OUT, Full_IN, Full_OUT, Mix_IN, Mix_OUT
REVEdgeFiltVel = zeros(6,length(RT))*NaN; 
WaveCoherence = zeros(2,length(RT))*NaN; 
BendDur = zeros(1,length(RT))*NaN;

for i = 1:length(RT)
    temp1 = Waves_Pointer{i,1};
    temp2 = Waves_Pointer{i,2};
    FRWVtemp = [];FRWVfullTemp = [];FRWVmixTemp = [];
    FRWVtempIN = [];
    FRWVtempOUT = [];
    REVVtemp = [];REVVfullTemp = [];REVVmixTemp = [];
    FRWfulltemp = []; REVfulltemp = [];
    BendDurTemp = [];
    FRWEdgeFiltT = cell(6,1);
    REVEdgeFiltT = cell(6,1);
    if size(temp1,1) ~= 0
        for j = 1:size(temp1,1)
            BendDurTemp = [BendDurTemp Waves1{temp1(j,1)}.Dur];
            if ~isempty(Waves1{temp1(j,1)}.FRWI)
                for k = 1:size(Waves1{temp1(j,1)}.FRWI,1)
                    if temp1(j,2) >= Waves1{temp1(j,1)}.FRWI(k,1) && temp1(j,2) <= Waves1{temp1(j,1)}.FRWI(k,2)
                        FRWVtemp = [FRWVtemp Waves1{temp1(j,1)}.FRWV(k)];
                        if Waves1{temp1(j,1)}.EdgePrc < 0.5
                            FRWEdgeFiltT{1} = [FRWEdgeFiltT{1} Waves1{temp1(j,1)}.FRWV(k)];
                        else
                            FRWEdgeFiltT{2} = [FRWEdgeFiltT{2} Waves1{temp1(j,1)}.FRWV(k)];
                        end
                        if Waves1{temp1(j,1)}.FRWF == 1
                            FRWfulltemp = [FRWfulltemp 1];
                            FRWVfullTemp = [FRWVfullTemp Waves1{temp1(j,1)}.FRWV(k)];
                            if Waves1{temp1(j,1)}.EdgePrc < 0.5
                                FRWEdgeFiltT{3} = [FRWEdgeFiltT{3} Waves1{temp1(j,1)}.FRWV(k)];
                            else
                                FRWEdgeFiltT{4} = [FRWEdgeFiltT{4} Waves1{temp1(j,1)}.FRWV(k)];
                            end
                        else
                            FRWfulltemp = [FRWfulltemp 0];
                            FRWVmixTemp = [FRWVmixTemp Waves1{temp1(j,1)}.FRWV(k)];
                            if Waves1{temp1(j,1)}.EdgePrc < 0.5
                                FRWEdgeFiltT{5} = [FRWEdgeFiltT{5} Waves1{temp1(j,1)}.FRWV(k)];
                            else
                                FRWEdgeFiltT{6} = [FRWEdgeFiltT{6} Waves1{temp1(j,1)}.FRWV(k)];
                            end
                        end
                    end
                end
            end
            if ~isempty(Waves1{temp1(j,1)}.REVI)
                for k = 1:size(Waves1{temp1(j,1)}.REVI,1)
                    if temp1(j,2) >= Waves1{temp1(j,1)}.REVI(k,1) && temp1(j,2) <= Waves1{temp1(j,1)}.REVI(k,2)
                        REVVtemp = [REVVtemp Waves1{temp1(j,1)}.REVV(k)];
                        if Waves1{temp1(j,1)}.EdgePrc < 0.5
                            REVEdgeFiltT{1} = [REVEdgeFiltT{1} Waves1{temp1(j,1)}.REVV(k)];
                        else
                            REVEdgeFiltT{2} = [REVEdgeFiltT{2} Waves1{temp1(j,1)}.REVV(k)];
                        end
                        if Waves1{temp1(j,1)}.REVF == 1
                            REVfulltemp = [REVfulltemp 1];
                            REVVfullTemp = [REVVfullTemp Waves1{temp1(j,1)}.REVV(k)];
                            if Waves1{temp1(j,1)}.EdgePrc < 0.5
                                REVEdgeFiltT{3} = [REVEdgeFiltT{3} Waves1{temp1(j,1)}.REVV(k)];
                            else
                                REVEdgeFiltT{4} = [REVEdgeFiltT{4} Waves1{temp1(j,1)}.REVV(k)];
                            end
                        else
                            REVfulltemp = [REVfulltemp 0];
                            REVVmixTemp = [REVVmixTemp Waves1{temp1(j,1)}.REVV(k)];
                            if Waves1{temp1(j,1)}.EdgePrc < 0.5
                                REVEdgeFiltT{5} = [REVEdgeFiltT{5} Waves1{temp1(j,1)}.REVV(k)];
                            else
                                REVEdgeFiltT{6} = [REVEdgeFiltT{6} Waves1{temp1(j,1)}.REVV(k)];
                            end
                        end
                    end
                end
            end
            
        end
    end
    if size(temp2,1) ~= 0
        for j = 1:size(temp2,1)
            BendDurTemp = [BendDurTemp Waves2{temp2(j,1)}.Dur];
            if ~isempty(Waves2{temp2(j,1)}.FRWI)
                for k = 1:size(Waves2{temp2(j,1)}.FRWI,1)
                    if temp2(j,2) >= Waves2{temp2(j,1)}.FRWI(k,1) && temp2(j,2) <= Waves2{temp2(j,1)}.FRWI(k,2)
                        FRWVtemp = [FRWVtemp Waves2{temp2(j,1)}.FRWV(k)];
                        if Waves2{temp2(j,1)}.EdgePrc < 0.5
                            FRWEdgeFiltT{1} = [FRWEdgeFiltT{1} Waves2{temp2(j,1)}.FRWV(k)];
                        else
                            FRWEdgeFiltT{2} = [FRWEdgeFiltT{2} Waves2{temp2(j,1)}.FRWV(k)];
                        end
                        if Waves2{temp2(j,1)}.FRWF == 1
                            FRWfulltemp = [FRWfulltemp 1];
                            FRWVfullTemp = [FRWVfullTemp Waves2{temp2(j,1)}.FRWV(k)];
                            if Waves2{temp2(j,1)}.EdgePrc < 0.5
                                FRWEdgeFiltT{3} = [FRWEdgeFiltT{3} Waves2{temp2(j,1)}.FRWV(k)];
                            else
                                FRWEdgeFiltT{4} = [FRWEdgeFiltT{4} Waves2{temp2(j,1)}.FRWV(k)];
                            end
                        else
                            FRWfulltemp = [FRWfulltemp 0];
                            FRWVmixTemp = [FRWVmixTemp Waves2{temp2(j,1)}.FRWV(k)];
                            if Waves2{temp2(j,1)}.EdgePrc < 0.5
                                FRWEdgeFiltT{5} = [FRWEdgeFiltT{5} Waves2{temp2(j,1)}.FRWV(k)];
                            else
                                FRWEdgeFiltT{6} = [FRWEdgeFiltT{6} Waves2{temp2(j,1)}.FRWV(k)];
                            end
                        end
                    end
                end
            end
            if ~isempty(Waves2{temp2(j,1)}.REVI)
                for k = 1:size(Waves2{temp2(j,1)}.REVI,1)
                    if temp2(j,2) >= Waves2{temp2(j,1)}.REVI(k,1) && temp2(j,2) <= Waves2{temp2(j,1)}.REVI(k,2)
                        REVVtemp = [REVVtemp Waves2{temp2(j,1)}.REVV(k)];
                        if Waves2{temp2(j,1)}.EdgePrc < 0.5
                            REVEdgeFiltT{1} = [REVEdgeFiltT{1} Waves2{temp2(j,1)}.REVV(k)];
                        else
                            REVEdgeFiltT{2} = [REVEdgeFiltT{2} Waves2{temp2(j,1)}.REVV(k)];
                        end
                        if Waves2{temp2(j,1)}.REVF == 1
                            REVfulltemp = [REVfulltemp 1];
                            REVVfullTemp = [REVVfullTemp Waves2{temp2(j,1)}.REVV(k)];
                            if Waves2{temp2(j,1)}.EdgePrc < 0.5
                                REVEdgeFiltT{3} = [REVEdgeFiltT{3} Waves2{temp2(j,1)}.REVV(k)];
                            else
                                REVEdgeFiltT{4} = [REVEdgeFiltT{4} Waves2{temp2(j,1)}.REVV(k)];
                            end
                        else
                            REVfulltemp = [REVfulltemp 0];
                            REVVmixTemp = [REVVmixTemp Waves2{temp2(j,1)}.REVV(k)];
                            if Waves2{temp2(j,1)}.EdgePrc < 0.5
                                REVEdgeFiltT{5} = [REVEdgeFiltT{5} Waves2{temp2(j,1)}.REVV(k)];
                            else
                                REVEdgeFiltT{6} = [REVEdgeFiltT{6} Waves2{temp2(j,1)}.REVV(k)];
                            end
                        end
                    end
                end
            end
        end
    end
    FRWVelocityRT(1,i) = mean(FRWVtemp);
    FRWVelocityRT(2,i) = mean(FRWVfullTemp);
    FRWVelocityRT(3,i) = mean(FRWVmixTemp);
    
    REVVelocityRT(1,i) = mean(REVVtemp);
    REVVelocityRT(2,i) = mean(REVVfullTemp);
    REVVelocityRT(3,i) = mean(REVVmixTemp);

    WaveCoherence(1,i) = sum(FRWfulltemp)/length(FRWfulltemp);
    WaveCoherence(2,i) = sum(REVfulltemp)/length(REVfulltemp);
    
    for j = 1:6
        FRWEdgeFiltVel(j,i) = mean(FRWEdgeFiltT{j});
        REVEdgeFiltVel(j,i) = mean(REVEdgeFiltT{j});
    end
    
    BendDur(i) = mean(BendDurTemp);
end

% Determine when bends and bends corresponding to behavior start for
% frequency analysis
BendStarts = zeros(1,length(RT));
FRWStarts = zeros(3,length(RT)); %All, Full, Mix
REVStarts = zeros(3,length(RT));
FRWStartsEdgeFilt = zeros(6,length(RT)); %All_IN, All_OUT, Full_IN, Full_OUT, Mix_IN, Mix_OUT
REVStartsEdgeFilt = zeros(6,length(RT));

for i = 1:length(Waves1)
    BendStarts(Waves1{i}.T(1)) = BendStarts(Waves1{i}.T(1)) + 1;
    if ~isempty(Waves1{i}.FRWI)
        for j = 1:length(Waves1{i}.FRWI(:,1))
            StartTime = Waves1{i}.FRWI(j,1);
            FRWStarts(1,Waves1{i}.T(StartTime)) = FRWStarts(1,Waves1{i}.T(StartTime)) + 1;
            if Waves1{i}.EdgePrc < 0.5
                FRWStartsEdgeFilt(1,Waves1{i}.T(StartTime)) = FRWStartsEdgeFilt(1,Waves1{i}.T(StartTime)) + 1;
            else
                FRWStartsEdgeFilt(2,Waves1{i}.T(StartTime)) = FRWStartsEdgeFilt(2,Waves1{i}.T(StartTime)) + 1;
            end
            
            if Waves1{i}.FRWF == 1
                FRWStarts(2,Waves1{i}.T(StartTime)) = FRWStarts(2,Waves1{i}.T(StartTime)) + 1;
                if Waves1{i}.EdgePrc < 0.5
                    FRWStartsEdgeFilt(3,Waves1{i}.T(StartTime)) = FRWStartsEdgeFilt(3,Waves1{i}.T(StartTime)) + 1;
                else
                    FRWStartsEdgeFilt(4,Waves1{i}.T(StartTime)) = FRWStartsEdgeFilt(4,Waves1{i}.T(StartTime)) + 1;
                end
            else
                FRWStarts(3,Waves1{i}.T(StartTime)) = FRWStarts(3,Waves1{i}.T(StartTime)) + 1;
                if Waves1{i}.EdgePrc < 0.5
                    FRWStartsEdgeFilt(5,Waves1{i}.T(StartTime)) = FRWStartsEdgeFilt(5,Waves1{i}.T(StartTime)) + 1;
                else
                    FRWStartsEdgeFilt(6,Waves1{i}.T(StartTime)) = FRWStartsEdgeFilt(6,Waves1{i}.T(StartTime)) + 1;
                end
            end
        end
    end
    if ~isempty(Waves1{i}.REVI)
        for j = 1:length(Waves1{i}.REVI(:,1))
            StartTime = Waves1{i}.REVI(j,1);
            REVStarts(1,Waves1{i}.T(StartTime)) = REVStarts(1,Waves1{i}.T(StartTime)) + 1;
            if Waves1{i}.EdgePrc < 0.5
                REVStartsEdgeFilt(1,Waves1{i}.T(StartTime)) = REVStartsEdgeFilt(1,Waves1{i}.T(StartTime)) + 1;
            else
                REVStartsEdgeFilt(2,Waves1{i}.T(StartTime)) = REVStartsEdgeFilt(2,Waves1{i}.T(StartTime)) + 1;
            end
            
            if Waves1{i}.REVF == 1
                REVStarts(2,Waves1{i}.T(StartTime)) = REVStarts(2,Waves1{i}.T(StartTime)) + 1;
                if Waves1{i}.EdgePrc < 0.5
                    REVStartsEdgeFilt(3,Waves1{i}.T(StartTime)) = REVStartsEdgeFilt(3,Waves1{i}.T(StartTime)) + 1;
                else
                    REVStartsEdgeFilt(4,Waves1{i}.T(StartTime)) = REVStartsEdgeFilt(4,Waves1{i}.T(StartTime)) + 1;
                end
            else
                REVStarts(3,Waves1{i}.T(StartTime)) = REVStarts(3,Waves1{i}.T(StartTime)) + 1;
                if Waves1{i}.EdgePrc < 0.5
                    REVStartsEdgeFilt(5,Waves1{i}.T(StartTime)) = REVStartsEdgeFilt(5,Waves1{i}.T(StartTime)) + 1;
                else
                    REVStartsEdgeFilt(6,Waves1{i}.T(StartTime)) = REVStartsEdgeFilt(6,Waves1{i}.T(StartTime)) + 1;
                end
            end
        end
    end
end
for i = 1:length(Waves2)
    BendStarts(Waves2{i}.T(1)) = BendStarts(Waves2{i}.T(1)) + 1;
    if ~isempty(Waves2{i}.FRWI)
        for j = 1:length(Waves2{i}.FRWI(:,1))
            StartTime = Waves2{i}.FRWI(j,1);
            FRWStarts(1,Waves2{i}.T(StartTime)) = FRWStarts(1,Waves2{i}.T(StartTime)) + 1;
            if Waves2{i}.EdgePrc < 0.5
                FRWStartsEdgeFilt(1,Waves2{i}.T(StartTime)) = FRWStartsEdgeFilt(1,Waves2{i}.T(StartTime)) + 1;
            else
                FRWStartsEdgeFilt(2,Waves2{i}.T(StartTime)) = FRWStartsEdgeFilt(2,Waves2{i}.T(StartTime)) + 1;
            end
            if Waves2{i}.FRWF == 1
                FRWStarts(2,Waves2{i}.T(StartTime)) = FRWStarts(2,Waves2{i}.T(StartTime)) + 1;
                if Waves2{i}.EdgePrc < 0.5
                    FRWStartsEdgeFilt(3,Waves2{i}.T(StartTime)) = FRWStartsEdgeFilt(3,Waves2{i}.T(StartTime)) + 1;
                else
                    FRWStartsEdgeFilt(4,Waves2{i}.T(StartTime)) = FRWStartsEdgeFilt(4,Waves2{i}.T(StartTime)) + 1;
                end
            else
                FRWStarts(3,Waves2{i}.T(StartTime)) = FRWStarts(3,Waves2{i}.T(StartTime)) + 1;
                if Waves2{i}.EdgePrc < 0.5
                    FRWStartsEdgeFilt(5,Waves2{i}.T(StartTime)) = FRWStartsEdgeFilt(5,Waves2{i}.T(StartTime)) + 1;
                else
                    FRWStartsEdgeFilt(6,Waves2{i}.T(StartTime)) = FRWStartsEdgeFilt(6,Waves2{i}.T(StartTime)) + 1;
                end
            end
        end
    end
    if ~isempty(Waves2{i}.REVI)
        for j = 1:length(Waves2{i}.REVI(:,1))
            StartTime = Waves2{i}.REVI(j,1);
            REVStarts(1,Waves2{i}.T(StartTime)) = REVStarts(1,Waves2{i}.T(StartTime)) + 1;
            if Waves2{i}.EdgePrc < 0.5
                REVStartsEdgeFilt(1,Waves2{i}.T(StartTime)) = REVStartsEdgeFilt(1,Waves2{i}.T(StartTime)) + 1;
            else
                REVStartsEdgeFilt(2,Waves2{i}.T(StartTime)) = REVStartsEdgeFilt(2,Waves2{i}.T(StartTime)) + 1;
            end
            if Waves2{i}.REVF == 1
                REVStarts(2,Waves2{i}.T(StartTime)) = REVStarts(2,Waves2{i}.T(StartTime)) + 1;
                if Waves2{i}.EdgePrc < 0.5
                    REVStartsEdgeFilt(3,Waves2{i}.T(StartTime)) = REVStartsEdgeFilt(3,Waves2{i}.T(StartTime)) + 1;
                else
                    REVStartsEdgeFilt(4,Waves2{i}.T(StartTime)) = REVStartsEdgeFilt(4,Waves2{i}.T(StartTime)) + 1;
                end
            else
                REVStarts(3,Waves2{i}.T(StartTime)) = REVStarts(3,Waves2{i}.T(StartTime)) + 1;
                if Waves2{i}.EdgePrc < 0.5
                    REVStartsEdgeFilt(5,Waves2{i}.T(StartTime)) = REVStartsEdgeFilt(5,Waves2{i}.T(StartTime)) + 1;
                else
                    REVStartsEdgeFilt(6,Waves2{i}.T(StartTime)) = REVStartsEdgeFilt(6,Waves2{i}.T(StartTime)) + 1;
                end
            end
        end
    end
end


clear('temp1','temp2','FRWVtemp','REVVtemp','FRWfulltemp','REVfulltemp','BendDurTemp','FRWVfullTemp','FRWVmixTemp','REVVfullTemp','REVVmixTemp','Range','StartTime');
clear('FRWEdgeFiltT','REVEdgeFiltT');

%%
save([Worm ' Behavior'],'SlidWindow','frameRate','scrsz','Worm','SegNum','NumFrames','SEGMENTS','frames','AHA_S05','RT','RT_L','IN','OUT',...
    'QBraw_RT','Qbouts_L','Mbouts_L','QboutsFilt','Mbouts','QB_RT','QB_RT_SumSeg',...
    'QuiescAll_LRT','MotAll_L','QboutAll','QboutAll_LRT','QboutAll_L','MboutAll_LRT',...
    'MboutAll_L','MboutAll','QBoutData','QAllBoutData','MBoutData','MAllBoutData',...
    'QBnextMB','MBnextQB','QBoutDur','QBoutNum','QCollective','CumSumSegQuiesc','CumSumNumSegQuiesc',...
    'Waves1','Waves2','Waves_Pointer','F','R','D','Q','Unc','Behavior',...
    'FRWVelocityRT','REVVelocityRT','FRWEdgeFiltVel','REVEdgeFiltVel','WaveCoherence','BendDur',...
    'BendStarts','FRWStarts','REVStarts','FRWStartsEdgeFilt','REVStartsEdgeFilt',...
    'EdgePrc','WormLengths');

save([Worm ' QuiescenceData'],'Worm','RT_L','QBraw_RT','Qbouts_L','Mbouts_L','QboutsFilt','Mbouts','QB_RT',...
    'QuiescAll_LRT','MotAll_L','QboutAll','QboutAll_LRT','QboutAll_L','MboutAll_LRT',...
    'MboutAll_L','MboutAll','QCollective','CumSumSegQuiesc','CumSumNumSegQuiesc');

save([Worm ' BendsData'],'Worm','RT_L','Waves1','Waves2','Waves_Pointer','F','R','D','Q','Unc','Behavior');





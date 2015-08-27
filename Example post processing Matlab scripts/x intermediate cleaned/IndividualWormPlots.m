
%%
cJet = colormap(jet(SegNum));
c1b = colormap(winter(SegNum));
c1c = colormap(autumn(SegNum));
close all;
scrsz = get(0,'screensize');
HourXaxis = RT/10/60/60;
getsize = @(x)x(1);
%% Plot Angles
figure;hold on;
for i = SEGMENTS
    plot((frames(FRAMES)),AHA_S05{5}(i,FRAMES),'Color',cJet(i,:))
    hold on;
end
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' Angles.fig']);
close;

%% Qbouts
%Plot the smoothed and identified quiescence
figure; hold on;
for j = SEGMENTS
    plot(frames(FRAMES),AHA_S05{5}(j,FRAMES),'Color',cJet(j,:))
    for i = 1:getsize(size(QboutAll{1}))
        plot(frames(QboutAll{1}(i,1):QboutAll{1}(i,(end))),AHA_S05{5}(j,QboutAll{1}(i,1):QboutAll{1}(i,(end))),'.','MarkerSize',10,'Color',c1b(j,:));
    end
end
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' AnglesWithAllQ.fig']);
close;

%% Bends
figure;hold on;
for i = 1:length(Waves1)
   plot(Waves1{i}.T,Waves1{i}.S,'k.','MarkerSize',10); plot(Waves1{i}.T(1),Waves1{i}.S(1),'k.','MarkerSize',30);
   plot(Waves1{i}.FRWT,Waves1{i}.FRWS,'g.','MarkerSize',40);
   plot(Waves1{i}.REVT,Waves1{i}.REVS,'r.','MarkerSize',30);
   plot(Waves1{i}.QuiescT,Waves1{i}.QuiescS,'c.','MarkerSize',20);
end
for i = 1:length(Waves2)
    plot(Waves2{i}.T,Waves2{i}.S,'k.','MarkerSize',10); plot(Waves2{i}.T(1),Waves2{i}.S(1),'k.','MarkerSize',30);
    plot(Waves2{i}.FRWT,Waves2{i}.FRWS,'g.','MarkerSize',40);
    plot(Waves2{i}.REVT,Waves2{i}.REVS,'r.','MarkerSize',30);
    plot(Waves2{i}.QuiescT,Waves2{i}.QuiescS,'c.','MarkerSize',20);
end
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' Waves.fig']);
close

%% Quiescence Analysis
figure;hold on;
plot(HourXaxis,QFSW(1,:),'k')
plot(HourXaxis,QFSW(2,:),'b')
ylim([0 1]);
xlabel('Time (hr)')
ylabel('Quiescence Fraction')
title([Worm ' Quiescence Fraction']);
legend('Standard','Body')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' QF.fig']);
saveas(gcf,[Worm ' QF.pdf']);
close;

figure;hold on;
plot(HourXaxis,QBoutDurSW(1,:),'k')
plot(HourXaxis,QBoutDurSW(2,:),'b')
xlabel('Time (hr)')
ylabel('Qbout Duration (sec)')
title([Worm ' Qbout Duration']);
legend('Standard','Body')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' QboutDur.fig']);
saveas(gcf,[Worm ' QboutDur.pdf']);
close;

%%
figure;hold on;
plot(HourXaxis,BehaviorPrc10min(1,:),'g');
plot(HourXaxis,BehaviorPrc10min(2,:),'r');
plot(HourXaxis,BehaviorPrc10min(3,:),'k');
plot(HourXaxis,BehaviorPrc10min(4,:),'b');
ylim([0 1]);
xlabel('Time (hr)')
ylabel('Fraction')
title([Worm ' Behavior Percentage']);
legend('FRW','REV','DWELL','QUIESC')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BehaviorPrc10min.fig']);
saveas(gcf,[Worm ' BehaviorPrc10min.pdf']);
close;

figure;hold on;
plot(HourXaxis,FRWVelocities(1,:),'k');
plot(HourXaxis,FRWVelocities(2,:),'r');
plot(HourXaxis,FRWVelocities(3,:),'b');
xlabel('Time (hr)')
ylabel('Velocity (%worm/sec)')
title([Worm ' Forward Velocities']);
legend('All','Full','Mix')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' FRWvelocities.fig']);
saveas(gcf,[Worm ' FRWvelocities.pdf']);
close;

figure;hold on;
plot(HourXaxis,REVVelocities(1,:),'k');
plot(HourXaxis,REVVelocities(2,:),'r');
plot(HourXaxis,REVVelocities(3,:),'b');
xlabel('Time (hr)')
ylabel('Velocity (%worm/sec)')
title([Worm ' Reverse Velocities']);
legend('All','Full','Mix')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' REVvelocities.fig']);
saveas(gcf,[Worm ' REVvelocities.pdf']);
close;

figure;hold on;
plot(HourXaxis,FRWFrequencies(1,:),'k');
plot(HourXaxis,FRWFrequencies(2,:),'r');
plot(HourXaxis,FRWFrequencies(3,:),'b');
xlabel('Time (hr)')
ylabel('Frequency (1/sec)')
title([Worm ' Forward Wave Frequencies']);
legend('All','Full','Mix')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' FRWfrequencies.fig']);
saveas(gcf,[Worm ' FRWfrequencies.pdf']);
close;

figure;hold on;
plot(HourXaxis,REVFrequencies(1,:),'k');
plot(HourXaxis,REVFrequencies(2,:),'r');
plot(HourXaxis,REVFrequencies(3,:),'b');
xlabel('Time (hr)')
ylabel('Frequency (1/sec)')
title([Worm ' Reverse Wave Frequencies']);
legend('All','Full','Mix')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' REVfrequencies.fig']);
saveas(gcf,[Worm ' REVfrequencies.pdf']);
close;

figure;
plot(HourXaxis,BendInitFrequency,'k');
xlabel('Time (hr)')
ylabel('Frequency (1/sec)')
title([Worm ' Bend Initiation Frequencies']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BendFrequencies.fig']);
saveas(gcf,[Worm ' BendFrequencies.pdf']);
close;

figure;
plot(HourXaxis,BendDuration,'k');
xlabel('Time (hr)')
ylabel('Duration (sec)')
title([Worm ' Bend Duration']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BendDuration.fig']);
saveas(gcf,[Worm ' BendDuration.pdf']);
close;

%%

figure;hold on;
plot(HourXaxis,WaveCoherenceSW(1,:),'k');
plot(HourXaxis,WaveCoherenceSW(2,:),'r');
xlabel('Time (hr)')
ylabel('Fraction of Waves Coherent')
title([Worm ' Bend Propagation Coherence']);
legend('FRW','REV')
set(gcf,'position',scrsz);
ylim([-0.2 1.2])
saveas(gcf,[Worm ' WaveCoherence.fig']);
saveas(gcf,[Worm ' WaveCoherence.pdf']);
close;

%%
% figure;hold on;
% plot(frames(FRAMES)/10/60/60,WormLengths(FRAMES),'k')
% plot(HourXaxis,WL_Slid,'r');
% xlabel('Time (hr)')
% ylabel('Length (pixels)')
% title([Worm ' Worm Length']);
% set(gcf,'position',scrsz);
% saveas(gcf,[Worm ' WormLength.fig']);
% saveas(gcf,[Worm ' WormLength.pdf']);
% close;


%%
% figure;hold on;
% plot(HourXaxis,AngSumSlid,'k');
% xlabel('Time (hr)')
% ylabel('Angle Sum (Radians)')
% title([Worm ' Angle Sum']);
% set(gcf,'position',scrsz);
% saveas(gcf,[Worm ' AngSum.fig']);
% saveas(gcf,[Worm ' AngSum.pdf']);
% close;

figure;hold on;
plot(HourXaxis,AngSumAbsSlid,'k');
xlabel('Time (hr)')
ylabel('abs(Angle) Sum (Radians)')
title([Worm ' abs(Angle) Sum']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' AngSumAbs.fig']);
saveas(gcf,[Worm ' AngSumAbs.pdf']);
close;

%%
figure; hold on;
for i = SEGMENTS
plot(HourXaxis,DAng(i,:),'Color',cJet(i,:))
end
plot(HourXaxis,MeanDAng,'k','LineWidth',2)
xlabel('Time (hr)')
ylabel('diff(Angle)')
set(gcf,'position',scrsz);
title([Worm ' DAng']);
saveas(gcf,[Worm ' DAng.fig']);
saveas(gcf,[Worm ' DAng.pdf']);
close;


figure;hold on;
plot(HourXaxis,CoM_AreaSW,'k');
xlabel('Time (hr)')
ylabel('Center of Mass area (pixels)')
title([Worm ' Center of Mass motion in 60 sec']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' CoM.fig']);
saveas(gcf,[Worm ' CoM.pdf']);
close;
















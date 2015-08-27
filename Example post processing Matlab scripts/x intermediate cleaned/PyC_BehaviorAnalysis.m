try
    cd('D:\AB_PyCelegans\N2_150313')
    load('N2_150313.mat');
    FTF = [1	137 %Frames to flip
        8025	14432
        14465	16970
        17159	29170
        41295	41849
        53880	54554
        54591	54710
        64888	65059
        65531	68512];
    FTD = [68513	68609]; %frames to delete
    Worm = 'N2_150313';
    NPath = 'D:\AB_PyCelegans\N2_150313\Analysis';
    mkdir(NPath);
    cd(NPath)
    NM = 0;
    N = 1;
    frameRate = 5;
    CalculateAngles_edit;
    BehaviorID_edit;
    CalculateSlidWindow_edit;
    IndividualWormPlots_edit;
catch err
    cd
    err.message
end
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUN resave.m first %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    cd('W:\PyCelegans Data\AgarExperiments\W70241_egl3')
    load('egl3_70241_N.mat');
    %%%%%%%%%%%% Head tail flips %%%%%%%%%%%%%%%%%%%%%%%%%%%
    FTF = [ % copy manualy head/tail flip indicies (two columns)
        
    ];  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%% Delete end frames to terminate exp early %%%%%%%%
    FTD = [ % manually enter to-delete-from frame no.
            % and LAST mmovie frame no. 
            9582, 12617; 
            22050, 27743; 
            31143, 34082
        ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Worm = 'egl3 70241 N';
    NPath = 'W:\PyCelegans Data\AgarExperiments\W70241_egl3\NAnalysis';
    mkdir(NPath);
    cd(NPath)
    % NM = 0; % not used 
    N = 1;
    frameRate = 5;
    CalculateAngles_edit;
    BehaviorID_edit;
    CalculateSlidWindow_edit;
    IndividualWormPlots_edit;
    EpochCalculations;
catch err
    cd
    err.message
end
clear
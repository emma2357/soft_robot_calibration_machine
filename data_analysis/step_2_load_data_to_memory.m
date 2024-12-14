% This script is used to load the results of multiple tests into a single
% set of arrays. These arrays are loaded onto the workspace and can be used
% by the step_3 script.


%-----------------------------------------
% Step 1 - Select the files containing the condensed testing data.
%          These files should be created using the "prepare_data" 
%          script.  
%-----------------------------------------

% create a small dataset for quick testing that contains only a fraction of
% the total data. This is an additionional data set (small_...) beyond 
% the original entire data set.
percentOfDataToUse = 0.2;   

% List all of the experiments you would like included.
experimentDirs      = [
    %"2024-10-09_12_07_23"; % curvature 10
    %"./Results/Force/", "2024-10-06_14_48_07"; -- old data
    %"2024-10-28_12_24_15"; % curvature 0, 5
    %"2024-10-30_11_57_53"; % all bending, most recent
    %"2024-10-08_14_50_05"; % all bending, old
    %"2024-09-21_14_45_41"; %-- old data 
    %"2024-11-13_14_57_22";  % new mismatched force data
    "2024-11-19_09_41_33"; % newest curvature data
    "2024-11-20_11_00_00"; % new force data 1
    "2024-11-21_12_21_36"; % new force data 2
    ];

%% load data
forceArray = [];
curvaOrienSensorArray = [];
sensorArray = [];
curvArray = [];
orienArray = [];
orienOneHotArray = [];
allInputArray = [];
allDataTable = [];

dsXTrain = [];
dsT2Train = [];
dsXTest = [];
dsT2Test = [];
dsXVal = [];
dsT2Val = [];

for fileInd = 1:length(experimentDirs)
    file_dir  = strcat("./data/", experimentDirs(fileInd, 1), ".mat");
    load(file_dir)

    forceArray_current = processed_data.orig.forces;
    forceArray         = [forceArray; forceArray_current];

    sensorArray_current = processed_data.orig.sensor_data;
    sensorArray         = [sensorArray; sensorArray_current];

    curvaOrienSensorArray_current = processed_data.orig.all_input_data;
    curvaOrienSensorArray = [curvaOrienSensorArray; curvaOrienSensorArray_current];

    curvArray_current = [curvaOrienSensorArray_current(:, 1)];
    curvArray         = [curvArray; curvArray_current];   

    orienArray_current = [curvaOrienSensorArray_current(:, 2)];
    orienArray         = [orienArray; orienArray_current];
    
    orienOneHotArray_current = processed_data.classification.orientation_one_hot;
    orienOneHotArray         = [orienOneHotArray; orienOneHotArray_current];
    
    allDataArray_current = processed_data.orig.all_data;
    allDataTable = [allDataTable; allDataArray_current];

    idxTrain = processed_data.orig.idxTrain;
    idxTest  = processed_data.orig.idxTest;
    idxVal   = processed_data.orig.idxVal;

    dsXTrain  = [dsXTrain; curvaOrienSensorArray_current(idxTrain, :)];
    dsT2Train = [dsT2Train; forceArray_current(idxTrain, :)];

    dsXTest  = [dsXTest; curvaOrienSensorArray_current(idxTest, :)];
    dsT2Test = [dsT2Test; forceArray_current(idxTest, :)];

    dsXVal   = [dsXVal; curvaOrienSensorArray_current(idxVal, :)];
    dsT2Val  = [dsT2Val; forceArray_current(idxVal, :)];
end


%-----------------------------------------
% Step 2 - Orientation predictions are difficult at low curvature
%         (all orientations are the same for a straight instrument)
%          Here the orientation labels can be removed (all made to be zero)
%          for low curvature.
%-----------------------------------------
zeroCuvInd        = curvArray < 3;
zeroed_orienArray = orienArray(:);
zeroed_orienArray(zeroCuvInd) = 0;

%-----------------------------------------
% Step 3 - Randomly select a percentage of the data for use with regression
%          learner. When the amount of data is too large, computation 
%          times are intractable.
%-----------------------------------------

selectIndex               = rand(length(sensorArray),1) < percentOfDataToUse;
small_curv_orien_sensor_array        = curvaOrienSensorArray(selectIndex,:);
small_sensor_array        = sensorArray(selectIndex,:);
small_curvat_array        = curvArray(selectIndex,:);
small_orient_array        = orienArray(selectIndex,:);
small_zeroed_orient_array = zeroed_orienArray(selectIndex,:);
small_force_array         = forceArray(selectIndex,:);


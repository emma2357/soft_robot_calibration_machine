% -----------------------------
% This script is used to gather the results from multiple 
% orientations and curvatures within a folder and create
% a single consistent .mat file with all of the data.
%
% This file only needs to be run once. 
% -----------------------------

clear

% -----------------------------
% Step 1 - Select all of the folders to consolidate. A single .mat 
%          file will be created for each of the folders. The folder
%          name is entered in two parts. The first part 
%          includes up to the Force or Bending directory. The second
%          entry on the row indicates the particular test folder to be 
%          analyzed. 
% -----------------------------
experimentDirs      = [
    %"./Results/Force/","2024-10-09_12_07_23"; % curvature 10
    %"./Results/Force/", "2024-10-06_14_48_07"; %-- old data
    %"./Results/Force/", "2024-10-28_12_24_15"; % curvature 0, 5
    %"./Results/Bending/", "2024-10-30_11_57_53"; % all bending, most recent
    %"./Results/Bending/", "2024-10-08_14_50_05"; % all bending, old
    %"./Results/Bending/", "2024-09-21_14_45_41" %-- old data 
    %"./Results/Force/", "2024-11-13_14_57_22"; % new mismatched force data
    "./Results/Bending/", "2024-11-19_09_41_33"; % newest curvature data
    "./Results/Force/", "2024-11-20_11_00_00"; % newest force data 1
    "./Results/Force/", "2024-11-21_12_21_36"; % newest force data 2
    ];

% training / validation / testing data split 
train_val_test = [0.8 0.1 0.1];

% -----------------------------
% Step 2 - Data is split into structs
% -----------------------------
for fileInd = 1:length(experimentDirs(:,1))
    experimentRootDir  = strcat(experimentDirs(fileInd, 1), experimentDirs(fileInd, 2));

    experiments = {};
    
    used_load_types = [-1, 1];
    processed_data = struct();
    
    % Gather all of the data from all orientations
    % load all .mat files found in sub directories using the 
    % filename to determine the orientation
    tmp_list_bending_directories = dir(strcat(experimentRootDir, "/*/*.mat"));
    list_bending_directories =[];
    for i = 1:length(tmp_list_bending_directories)
        if (~startsWith(tmp_list_bending_directories(i).name, "."))
            list_bending_directories = [list_bending_directories; tmp_list_bending_directories(i)];
        end
    end
    
    nrExperiments = length(list_bending_directories);
    for k = 1:nrExperiments
        file_dir = list_bending_directories(k).folder;
        file_name = list_bending_directories(k).name;
        load(strcat(file_dir, "/", file_name));
        orientation = extractBefore(file_name, "_");
    
        experiments{end + 1} = struct("orientation", ...
            str2num(orientation), "fileName", file_name, ...
            "voltages", loss_table, "loss_table", loss_table, ...
            "extra_data", d);
    end
    
    % pull the sensor data and force data from each experiment
    % and concaternate into a single array.
    sensorReadingArray = [];
    voltageReadingArray= [];
    forceArray = [];
    curvaArray = [];
    orienArray = [];
    allDataArray = [];
    for i = 1:nrExperiments
        rawVoltageTable = experiments{i}.voltages;
        rawVoltageTable = rawVoltageTable(find(~isnan(rawVoltageTable.("F[N]"))), :);
    
        if exist('rawVoltageTable.loadType','var')
            voltageTable = rawVoltageTable(ismember(rawVoltageTable.loadType, used_load_types), :);
        else
            voltageTable = rawVoltageTable;
        end
    
        sensorReadingArray = [sensorReadingArray; voltageTable.loss0, voltageTable.loss1, voltageTable.loss2, voltageTable.loss3, voltageTable.loss4];
        voltageReadingArray = [voltageReadingArray; voltageTable.("A0[V]"), voltageTable.("A1[V]"), voltageTable.("A2[V]"), voltageTable.("A3[V]"), voltageTable.("A4[V]")];
        
        forceArray = [forceArray; voltageTable.("F[N]")];
        curvaArray = [curvaArray; voltageTable.("curvature[m^-1]")];
        for j = 1:numel(voltageTable.("F[N]"))
            orienArray = [orienArray;  experiments{i}.orientation];
        end
        
        allDataArray = [allDataArray; experiments{i}.loss_table];
        allDataArray = allDataArray(find(~isnan(allDataArray.("F[N]"))), :);
    end
    

    forceArray(forceArray < 0) = 0;
    forceArray(forceArray > 7) = 7;
    
    % merge datasets
    curvaOrienSensorArray = [curvaArray, orienArray, sensorReadingArray];
    allDataArray = [allDataArray, array2table(orienArray,"VariableNames","orientation")];

    totRev_X_Array = [forceArray, curvaArray, orienArray];
    totRev_Y_Array = [voltageReadingArray];

    % Partition used to split into training and test data
    n = length(sensorReadingArray);
    [idxTrain,idxTest, idxVal] = trainingPartitions(n,train_val_test);

    % save the original data
    orig = struct();
    orig.totRev_X_Array = totRev_X_Array(:,:);
    orig.totRev_Y_Array = totRev_Y_Array(:,:);
    orig.forces      = forceArray(:);
    orig.curvature   = curvaArray(:);
    orig.orientation = orienArray(:);
    orig.sensor_data = sensorReadingArray(:,:);
    orig.all_input_data = curvaOrienSensorArray(:,:);
    orig.all_data = allDataArray(:,:);
    orig.idxTrain = idxTrain;
    orig.idxTest = idxTest;
    orig.idxVal = idxVal;
    processed_data.orig = orig;
    
    %hpartition = cvpartition(n, 'Holdout', percent_test);
    %processed_data.hpartition = hpartition;
    
    % save (orientation) data for classification analysis
    classification = struct();
    
    % gives [1, 2, 3, 4, 5, 6] for [0, 60, 120, 180, 240, 300]
    classificationOrienArray = round(orienArray/60) + 1;
    
    classification.orientation         = classificationOrienArray(:);
    [classification.num_classes, classification.class_weights ] = computeClassWeights(classification.orientation);
    
    orienOneHotArray = oneHotEncodeArray(classificationOrienArray, classification.num_classes);

    classification.orientation_one_hot = orienOneHotArray;
    [classification.dsTrain, classification.dsVal] = splitIntoValidationAndTrainingData(orig.sensor_data, classification.orientation_one_hot, idxTrain,idxTest, idxVal);

    processed_data.classification = classification;


    % save data for regression analysis
    regression = struct();
    % Split into training testing and validation sets
    [regression.dsTrain, regression.dsTest, regression.dsVal] = ...
        splitIntoValidationAndTrainingData(orig.all_input_data, ...
        orig.forces, idxTrain,idxTest, idxVal);
    
    processed_data.regression = regression; 
    
    %% Save calibration data
    save(strcat("./data/", experimentDirs(fileInd, 2)), "processed_data");
end

% reorganize the arrays and create a datastore for training, testing, and validation.
function [dsTrain, dsTest, dsVal] = splitIntoValidationAndTrainingData(xData, yData, idxTrain,idxTest, idxVal)
    xTrain = xData(idxTrain, :);
    yTrain = yData(idxTrain, :);
    
    xTest = xData(idxTest, :);
    yTest = yData(idxTest, :);

    xVal = xData(idxVal, :);
    yVal = yData(idxVal, :);
    
    %% Create datastores out of the .mat file data
    % train
    dsXTrain = arrayDatastore(xTrain);
    dsT2Train = arrayDatastore(yTrain);
    
    % test
    dsXTest = arrayDatastore(xTest);
    dsT2Test = arrayDatastore(yTest);

    % validate
    dsXVal = arrayDatastore(xVal);
    dsT2Val = arrayDatastore(yVal);
    
    dsTrain = combine(dsXTrain, dsT2Train);
    dsTest = combine(dsXTest, dsT2Test);
    dsVal = combine(dsXVal, dsT2Val);
end

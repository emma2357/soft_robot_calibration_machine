% This file can be used to show the results 
% of a trained machine learning model applied to 
% data that is already loaded into memory.
% Use the "loadDataToMemory" script to gather 
% and filter data before executing this script. 

%-----------------------------------------
% STEP 1: Select plot parameters
%-----------------------------------------
largeFont     = 20;
smallFont     = 10;
dotSize       = 5;


%-----------------------------------------
% Step 2: Load a model used to predict force, curvature and orienation.
%         
%-----------------------------------------
use7InputForForce = true;
showCurv  = true;
showOrien = true;
showForce = true;
predictOrienCurvforForce = true;

% 5 input force models
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_5inputs_nn_32_32_tanh.mat").trainedModel;
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_5inputs_gau_ratQuad.mat").trainedModel;
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_5inputs_gau_expGPR.mat").trainedModel;
%force_nn = load("./TrainedModels/RL_20perdata_force_7inputs_nn_32_32_relu.mat").trainedModel;

% 7 input force models
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_7inputs_nn_32_32_tanh.mat").trainedModel;
force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_7inputs_gau_ratQuad.mat").trainedModel;
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_7inputs_gau_expGPR.mat").trainedModel;

%curvature models
%curv_nn = load("./TrainedModels_20perdata/RL_20perdata_curv_nn_32_32_tanh.mat").trainedModel;
curv_nn = load("./TrainedModels_20perdata/RL_20perdata_curv_gau_expGPR.mat").trainedModel;
%curv_nn = load("./TrainedModels_20perdata/CL_20perdata_curv_knn_fine.mat").trainedModel;

% orientation models
%orien_nn = load("./TrainedModels_20perdata/RL_20perdata_orien_nn_32_32_relu.mat").trainedModel;
orien_nn = load("./TrainedModels_20perdata/CL_20perdata_orien_knn_fine.mat").trainedModel;

% Gaussian Regression
% 7 input force models
%force_nn = load(".\TrainedModels\regression_simple_nn.mat").net;
%force_nn = load(".\TrainedModels\regressionLearner_nn_7inputs_allCurvature_Bending_32_32.mat").trainedModel_nn;

%-----------------------------------------
% Step 3: Copy input data and replace curvature and orientation columns in pred data array
%         with predicted values.
%-----------------------------------------
dsXTrain_pred = dsXTrain(:,:);

if predictOrienCurvforForce
    if (hasPredictFcn(curv_nn))
        dsXTrain_pred(:,1) = curv_nn.predictFcn(dsXTrain(:,3:7));
    else
        dsXTrain_pred(:,1) = curv_nn.predict(dsXTrain(:,3:7));
    end
    
    if (hasPredictFcn(orien_nn))
        dsXTrain_pred(:,2) = orien_nn.predictFcn(dsXTrain(:,3:7));
    else
        dsXTrain_pred(:,2) = orien_nn.predict(dsXTrain(:,3:7));
    end
end

%-----------------------------------------
% Step 4: Plot the accuracy of the curvature  predictions.
%-----------------------------------------
clf;
if showCurv
    tiledlayout(2, 3);
    nexttile
    hold on;
    
    xlabel("True Curvature (1/m)", "FontSize", largeFont);
    ylabel("Predicted Curvature (1/m)", "FontSize", largeFont);
    title("Accuracy of Curvature Predictions", "FontSize", largeFont);
    grid on;
    plot(dsXTrain(:,1), dsXTrain_pred(:,1), '.', 'MarkerSize',dotSize);
    
    addPerfectFitPlusPercent([0:11], 0.2);
end

%-----------------------------------------
% Step 5: Plot the accuracy of the orientaiton predictions.
%-----------------------------------------

if showOrien
    nexttile
    hold on;
    xlabel("True Orientation (Deg)", "FontSize", largeFont);
    ylabel("Predicted Orientation (Deg)", "FontSize", largeFont);
    title("Accuracy of Orientation Predictions", "FontSize", largeFont);
    grid on;
    plot(dsXTrain(:,2), dsXTrain_pred(:,2), '.', 'MarkerSize',dotSize);
    
    addPerfectFitPlusPercent([0:370], 0.2);
end 

%-----------------------------------------
% Step 6: predict force values.
%-----------------------------------------
if showForce
    if (use7InputForForce)
        tmpXData = dsXTrain_pred;
    else
        tmpXData = dsXTrain(:, 3:7);
    end
    
    if (hasPredictFcn(force_nn))
       [Y1] = force_nn.predictFcn(tmpXData);
    else
       [Y1] = predict(force_nn, tmpXData);
    end
    
    dataY = dsT2Train;


%-----------------------------------------
% Step 7: Plot the force accuracy.
%-----------------------------------------
    %figure;
    %clf;
    %tiledlayout(2, 2);
    nexttile
    hold on
    fontsize(smallFont, "points");
    xlabel("True Force (N)", "FontSize", largeFont)
    ylabel("Predicted Force (N)", "FontSize", largeFont)
    title("Accuracy of Network's Predictions", "FontSize", largeFont)
    grid on
    plot(dataY, Y1, '.', 'MarkerSize',dotSize)
    
    addPerfectFitPlusPercent([0:7], 0.2);
    
    %print(gcf, 'predicted_vs_actual_force.png', '-dpng', '-r300')
    
%-----------------------------------------
% Step 8: split forces into high (5N+ and low 0-5N) and plot 
%         confusion matrix.
%-----------------------------------------
    f_newtons = sum(abs(dataY - Y1) < 0.5) / length(Y1)
    
    predicted_twoCategory = Y1 >= 5;
    actual_twoCategory    = dataY >= 5;
    
    nexttile;
    
    hold off 
    label = {'0-5N','5N+'}; 
    label = categorical(label)
    C_two_cat = confusionmat(actual_twoCategory,predicted_twoCategory)
    if (numel(C_two_cat) == 1)
        label = {'0-5N'};
    end
    
    
    cm_two_cat = confusionchart(C_two_cat, label, 'RowSummary','row-normalized','ColumnSummary','column-normalized')
    
    xlabel(cm_two_cat, 'Predicted');
    ylabel(cm_two_cat, 'True');
    title(cm_two_cat, "Two Category Confusion Matrix ");
    
    fontsize(largeFont, "points");
    
    %print(gcf, 'confusionMatrix2x2.png', '-dpng', '-r300')
    
%-----------------------------------------
% Step 9: split forces into high 0-3, 3-5, and 5+ and plot 
%         confusion matrix.
%-----------------------------------------
    nexttile
    hold off
    
    classify = @(x) (x < 3) * 1 + (x >= 3 & x < 5) * 2 + (x >= 5) * 3;
    
    % Classify predicted and actual data
    predicted_classes = arrayfun(classify, Y1);
    actual_classes = arrayfun(classify, dataY);
    
    % Initialize confusion matrix
    confusion_matrix = zeros(3, 3);
    
    % Populate confusion matrix
    for i = 1:length(predicted_classes)
        confusion_matrix(actual_classes(i), predicted_classes(i)) = confusion_matrix(actual_classes(i), predicted_classes(i)) + 1;
    end
    
    %disp(confusion_matrix)
    label = {'0-3','3-5', '5+'};
    label = categorical(label, label, 'Ordinal', true);
    cm_three_cat = confusionchart(confusion_matrix, label, 'RowSummary','row-normalized','ColumnSummary','column-normalized');
    
    fontsize(largeFont, "points");
    
    xlabel(cm_three_cat, 'Predicted Force (N)');
    ylabel(cm_three_cat, 'True Force (N)');
    title(cm_three_cat, 'Three Category Confusion Matrix');
end 


disp("Statistics: Two Category")
computeStats(C_two_cat)
disp("Statistics: Three Category")
%computeStats(confusion_matrix)

% plot bounding lines that are +/- a given percent from a value.
function addPerfectFitPlusPercent(xValues, percent)
    tp  = xValues + percent * xValues;
    bp  = xValues - percent * xValues;
    plot(xValues, xValues);
    plot(xValues, tp, 'b');
    plot(xValues, bp, 'b');
    xlim([min(xValues) max(xValues)]);
    ylim([min(xValues) max(xValues)]);
end

% hasPredictFcn - try to determine if the structure v has 
%                 the predictFcn member.
function b = hasPredictFcn(v)
    try
        a = v.predictFcn;
        b = true;
    catch
        b = false;
    end
end

% computeStats - compute true/false positive / negatives ffrom 
%                a confusion matrix.
function [tp, tn, fp, fn] = computeStats(CM)
    totNr = sum(CM, "all")
    T = trace(CM)
    disp(["Accuracy:", (T)/(totNr)]);
    for i=1:length(CM)
        TP(i)=CM(i,i);
        FP(i)=sum(CM(:,i))-CM(i,i);
        TN(i)=sum(CM(:))-sum(CM(i,:))-sum(CM(:,i))+CM(i,i);
        FN(i)=sum(CM(i,:))-CM(i,i);
        disp(['True  Positive (', num2str(i), '):', num2str(TP(i))])
        disp(['False Positive (', num2str(i), '):', num2str(FP(i))])
        disp(['True  Negative (', num2str(i), '):', num2str(TN(i))])
        disp(['False Negative (', num2str(i), '):', num2str(FN(i))])
    end
    
end
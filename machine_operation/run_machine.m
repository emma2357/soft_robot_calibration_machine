% -------------------------------------------------------------------
%  Code for running soft rotobic sleeve calibration machine. 
% -------------------------------------------------------------------

clear
clear global

global m u g d s

h = gcf;
set(h, "Position", [0 0 1000 1000])

% Create Save directories
mkdir("Results")
mkdir("Results\Bending")
mkdir("Results\Force")

% -------------------------------------------------------------------
% Step 1. Select User Settings
% u - USER SETTINGS
% -------------------------------------------------------------------
u = struct();

% Measure distance from rotation point to start of sleeve restraint
u.armLength                  = 0.095; % meters

% Select a list of curvatures between 0 and 10 (for force loading)
u.used_curvatures            = [10];

% 360 is maximum
u.used_bending_orientations  = [0, 60, 120, 180, 240, 300];
u.used_force_orientations    = [0, 60];

u.nr_load_cycles             = 3;
u.rest_time_between_cycles   = 2; % seconds 

u.used_angles          = 2.0 * atand(1.0 ./ (u.used_curvatures .* u.armLength));
u.used_motor_settings  = getArmPositionFromAngle(u.used_angles);

% this code needs a curvature & orientation model attached
% however, set this to false and the 7 input force model won't concatenate
% the curvature & orientation predictions!
u.predictCurvOrien = false;

% -------------------------------------------------------------------
% Step 2. Select trained networks to be used for making real time
%         curvature, force, and orientaiton predictions.
% -------------------------------------------------------------------
% 5 input force prediction models
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_5inputs_nn_32_32_tanh.mat").trainedModel;
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_5inputs_gau_ratQuad.mat").trainedModel;
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_5inputs_gau_expGPR.mat").trainedModel;
%force_nn = load("./TrainedModels/RL_20perdata_force_7inputs_nn_32_32_relu.mat").trainedModel;

% 7 input force prediction models
force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_7inputs_nn_32_32_tanh.mat").trainedModel;
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_7inputs_gau_ratQuad.mat").trainedModel;
%force_nn = load("./TrainedModels_20perdata/RL_20perdata_force_7inputs_gau_expGPR.mat").trainedModel;

%curvature prediction models
curv_nn = load("./TrainedModels_20perdata/RL_20perdata_curv_nn_32_32_tanh.mat").trainedModel;
%curv_nn = load("./TrainedModels_20perdata/RL_20perdata_curv_gau_expGPR.mat").trainedModel;
%curv_nn = load("./TrainedModels_20perdata/CL_20perdata_curv_knn_fine.mat").trainedModel;

% orientation prediction  models
%orien_nn = load("./TrainedModels_20perdata/RL_20perdata_orien_nn_32_32_relu.mat").trainedModel;
orien_nn = load("./TrainedModels_20perdata/CL_20perdata_orien_knn_fine.mat").trainedModel;

% ------------------------------------------------------------------
% g - Graphs
%     Create graphs that will be updated with real time data.
% ------------------------------------------------------------------
g = struct();
g.yAxisMin                 = 0;  % Newtons
g.yAxisMax                 = 10; % Newtons

tiledlayout(5, 1)
nexttile;
g.force_graph = plot([0], [0], '-r');
grid on;
xlabel("Time [s]");
ylabel("Force [N]");

nexttile;
g.loss_graph  = plot([0], [0], '-r', [0], [0], '-g', [0], [0], '-b', [0], [0], '-k', [0], [0], '-m');
grid on;
xlabel("Time [s]");
ylabel("Loss");

nexttile;
g.predictedForce_graph  = plot([0], [0], '-r', [0], [0], 'bx');
legend("True", "Predicted")
grid on;
xlabel("Time [s]");
ylabel("Force");

nexttile;
g.predictedCurv_graph  = plot([0], [0], '-r', [0], [0], 'bx');
legend("True", "Predicted")
grid on;
xlabel("Time [s]");
ylabel("Curvature");

nexttile;
g.predictedOrien_graph  = plot([0], [0], '-r', [0], [0], 'bx');
legend("True", "Predicted")
grid on;
xlabel("Time [s]");
ylabel("Orientation");

% ------------------------------------------------------------------
% d - Data Store
%     All data will be stored in this data structure during the 
%     machine runs.
% ------------------------------------------------------------------
d = struct();
d.experiments = {};
d.voltages = [];
d.baseline_values = [];
d.baselineArray = [];

% ------------------------------------------------------------------
% m - Machine Objects
% ------------------------------------------------------------------
% Pin Layout
% D2 / D3 - HX711 senor giving the force reading
% D5      - Motor Changing the Curvature PWM 
% D6      - Motor Moving the Load Cell PWM 
% D7      - Limit Switch Upper
% D4      - Limit Switch Lower

m = struct();
m.upper_limit_pin      = "D7";
m.lower_limit_pin      = "D4";
m.force_servo_pin      = "D6";
m.position_servo_pin   = "D5";
m.rotation_servo_pin_L = "D9";
m.rotation_servo_pin_R = "D10";

m.calibration_mass     = 200;   % calibration mass in grams
m.load_set_point_N     = 6.0;   % in Newtons - load system up to this point
m.contact_mass         = 50;    % reading used for determining sleeve contact 
                                
% Constants obtained from load cell calibration with a 200g mass.
m.tare_weight  = 8.386536146000000e+06;   
m.scale_factor = 7.980477550000046e+02;


% Used for running the code on the lab PC
m.arduino_obj     = arduino("com5", "Uno", "libraries", {"Servo", 'advancedHX711/advanced_HX711'});
% Used for running on a MAC 
%m.arduino_obj    = arduino("/dev/cu.usbmodem1101", "Uno", "libraries", {"Servo", 'advancedHX711/advanced_HX711'});

% 5 kg load cell signal read through a HX711 chip
m.load_cell       = addon(m.arduino_obj, "advancedHX711/advanced_HX711", ...
    "Pins", {"D2","D3"}, "Gain", 128);
m.calibration_obj = calibration(100, m.calibration_mass);
setCalibration(m.calibration_obj, m.tare_weight, m.scale_factor);

% Set the limit switch pins
configurePin(m.arduino_obj, m.upper_limit_pin, 'pullup');
configurePin(m.arduino_obj, m.lower_limit_pin, 'pullup');
 
% start the timer
tic

% GoBilda DC Motor controller 
m.force_servo_obj      = servo(m.arduino_obj, m.force_servo_pin,    ...
    'MinPulseDuration', 1050e-6, 'MaxPulseDuration', 1950e-6);

% GoBilda Torque 5-turn servo
m.position_servo_obj   = servo(m.arduino_obj, m.position_servo_pin,  ...
    'MinPulseDuration', 700e-6, 'MaxPulseDuration', 2300e-6);

% GoBilda Torque 1-turn servo
m.rotation_servo_obj_L = servo(m.arduino_obj, m.rotation_servo_pin_L, ...
    'MinPulseDuration', 500e-6, 'MaxPulseDuration', 2500e-6);
m.rotation_servo_obj_R = servo(m.arduino_obj, m.rotation_servo_pin_R, ...
    'MinPulseDuration', 500e-6, 'MaxPulseDuration', 2500e-6);

% using gear ratio to compute sleeve rotation degrees per degree rotation of servo
m.max_rotation_servo_angle = 300 * 60/48;  

% ------------------------------------------------------------------
% s - Machine settings (used to save current machine positions)
% ------------------------------------------------------------------
s = struct();

% change curvature to 0 (straight tube)
if exist('StoreSettings\machineSettings.mat', 'file') == 2
    load StoreSettings\machineSettings.mat s
else
    s.cur_curv_pos = 1.0;
    s.cur_orie_pos = 0.0;
end

setArmPosition(s.cur_curv_pos);
setFiberOrientation(s.cur_orie_pos * m.max_rotation_servo_angle);
pause(3.0)

moveForceToStartPosition()
resetVoltages()

while (true)
    disp('Beginning a new test. Press control-C to exit.');
    disp('Select an option:')
    disp('  F - Test Forces')
    disp('  C - Test Curvature')
    disp('  R - Reset system (move force / arm / fiber to start position)')
    disp('  M - Manual rotations')
    disp('  T - Move force system to contact')
    disp('  S - Move force system to start')
    disp('  V - Display current readings')
    %disp('  CAL - Calibrate load cell')
    disp('  Q - quit')
    
    option = input("Enter Option: ", "s");
    if (option == "C")
        dateTimeStr = string(datetime('now','Format','yyyy-MM-dd_HH_mm_ss'));
    
        for orientation = u.used_bending_orientations
            setFiberOrientation(orientation);
            pause(2);
            testCurvatures(dateTimeStr, orientation);
        end
    elseif (option == "F")
        dateTimeStr = string(datetime('now','Format','yyyy-MM-dd_HH_mm_ss'));
    
        % output current load cell reading
        initial_load = getLoadInNewtons(5);
        disp(['Current Load is: ' num2str(initial_load) ' [N]']);
        
        % check that load reading is reasonable (should be zero). 
        if (abs(initial_load) > 0.3)
            custom_error("Initial load is too high. Check the machine setup.")
        end

        for j = 1:length(u.used_motor_settings)
            disp("Begin testing all fiber orientations for the new curvature setting.")
            disp(["Current curvature being tested is " + num2str(u.used_curvatures(j))])
            
            disp("USER ACTION REQUIRED ------------------------------------------------")
            disp("Replace 3D printed Force Application Device. Press enter to continue.")
            pause()
            setArmPosition(1);
            setFiberOrientation(0);
            for orientation = u.used_force_orientations
                setFiberOrientation(orientation);
                pause(2);
                testForces(dateTimeStr, m.load_set_point_N, orientation, j);
            end
        end
    elseif (option == "R")
        moveForceToStartPosition()
        setArmPosition(1);
        setFiberOrientation(0);
    elseif (option == "M")
        inputArmAngle = -1;
        while (inputArmAngle > 180 || inputArmAngle < 90)
            inputArmAngle = input("Enter arm position (180 to 90 degrees):")
            inputArmPos = getArmPositionFromAngle(inputArmAngle);
            disp(["Motor setting is:" num2str(inputArmPos)])
        end
        setArmPosition(inputArmPos);
        pause(1);

        inputFiberPos = -1;
        while (inputFiberPos > m.max_rotation_servo_angle || inputFiberPos < 0)
            inputFiberPos = input("Enter fiber orientation (0 to 360 degrees): ")
        end
        setFiberOrientation(inputFiberPos);
        pause(1);
    elseif (option == "Q")
        custom_error("Quit.")
    elseif (option == "T")
        moveToContact()
    elseif (option == "S")
        moveForceToStartPosition()
    elseif (option == "V")
        disp(["Force [N]: " num2str(getLoadInNewtons(5))]);
        [v1, v2, v3, v4, v5] = readLossVoltages();
        disp(["Light [V]: " num2str(v1) num2str(v2) num2str(v3) num2str(v4) num2str(v5)])  
    elseif (option == "CAL")
        m.calibration_obj = calibration(100, m.calibration_mass);
        disp('Press a key to tare')
        pause;
        tare(m.calibration_obj, load_cell)
        disp('Load weight and press key')
        scale(m.calibration_obj,load_cell);
        disp('Press a key to continue')
        pause;
    else
        disp("Incorrect selection. Try again.")
    end
end

% testForces - complete a load / unload set of cycles for a given curvature
%              and orientation.
function testForces(dateTimeStr, maxLoad, orientation, curvatureIndex)
    global u d
    c = onCleanup(@() finalizePositions());

    % Stop force motor from moving
    stopMotor();
    
    % Create new directory with timestamp
    mkdir("Results/Force/" + dateTimeStr + "/");
    dirName = "Results/Force/" + dateTimeStr + "/" + num2str(orientation);
    [status, msg, msgID] = mkdir(dirName);

    % Reset the machine to the starting position
    moveForceToStartPosition()
    setArmPosition(1);

    % Record baseline voltages
    recordBaselineVoltages();    
    saveBaselineData(dirName, num2str(orientation), 0);
    
    % Find the next arm curvature to be tested
    positionServoSetting = u.used_motor_settings(curvatureIndex);
    curvature            = u.used_curvatures(curvatureIndex);

    % Change sleeve curvature 
    setArmPosition(positionServoSetting); 
    pause(1.5)
    
    d.voltages = [];
    
    disp("Move until force applicator is in contact with the sleeve.")
    moveToContact()
 
    % Cycle load / unload cycles 
    disp("Load Cycles started.")
    for i = 1:u.nr_load_cycles
        disp(['Start Load Cycle: ' num2str(i)])
        disp('Resting.');
        collectDataFor(u.rest_time_between_cycles, i, curvature, orientation);
        disp('Load Increasing')
        gotoLoad(i, true, maxLoad,  curvature, orientation);
        stopMotor();
        disp('Resting.');
        collectDataFor(u.rest_time_between_cycles, i, curvature, orientation);
        disp('Load Decreasing')
        gotoLoad(i, false, 0, curvature, orientation);
        stopMotor();
    end

    stopMotor();

    % Save data to a file
    saveNewFile(dirName, num2str(orientation), curvature, "force");
    resetVoltages();
    
    % reset the force axis
    moveForceToStartPosition()
end

function collectDataFor(t, cycle_id, curvature_val, orientation)
    t0 = toc;
    while ((toc - t0) < t)
        pause(0.05);  
        recordForceData(false, cycle_id, 0, 1, 1, curvature_val, orientation);  
    end
end

function setArmPosition(curvature_pos)
    global m s
    for setting = s.cur_curv_pos:0.005*sign(curvature_pos-s.cur_curv_pos):curvature_pos
        writePosition(m.position_servo_obj, setting);
        pause(0.1);
    end
    
    s.cur_curv_pos = curvature_pos;
    writePosition(m.position_servo_obj, curvature_pos);
    save StoreSettings\machineSettings.mat s -mat;
end

function setFiberOrientation(orientationAngle)
    global m s
    orientationPos = orientationAngle/m.max_rotation_servo_angle;

    for setting = s.cur_orie_pos:0.005*sign(orientationPos-s.cur_orie_pos):orientationPos
        writePosition(m.rotation_servo_obj_L, (1.0 - setting));
        writePosition(m.rotation_servo_obj_R, (0.0 + setting));
        pause(0.1);
    end
    
    s.cur_orie_pos = orientationPos;
    writePosition(m.rotation_servo_obj_L, (1.0 - orientationPos));
    writePosition(m.rotation_servo_obj_R, (0.0 + orientationPos));
    save StoreSettings\machineSettings.mat s -mat;
end

% recordBaselineVoltages - record 50 baseline intensity values
function recordBaselineVoltages()
    global d

    disp("Reading Baseline Voltages")
    d.baselineArray = [];
    
    initial_force = getLoadInNewtons(5);

    for i=1:1:50
        [v1, v2, v3, v4, v5] = readLossVoltages();
        d.baselineArray = [d.baselineArray; v1, v2, v3, v4, v5, initial_force, 0, 0];
    end
    
    d.baseline_values = [];
    for i = 1:5
        d.baseline_values(i) = mean(d.baselineArray(:,i));
    end
end

% testCurvatures - record intensity as the curvature of the sleeve is
%                  changed
function testCurvatures(dateTimeStr, orientation)
    global m u
    global motorSettingFromAngle 
    c=onCleanup(@() finalizePositions());

    % Move load cell to top of machine and stop when limit switch
    % is pressed. 
    moveForceToStartPosition()
    
    mkdir("Results/Bending/" + dateTimeStr + "/");
    dirName = "Results/Bending/" + dateTimeStr + "/" + num2str(orientation);
    [status, msg, msgID] = mkdir(dirName);
   
    recordBaselineVoltages();    
    saveBaselineData(dirName, num2str(orientation), 0);
    
    nrCurvatureTrials = 1;

    for i = 1:1:nrCurvatureTrials
        for angle_val = 180:-10:90
            if (angle_val < 0.5)
                curvature_val = 0;
            else
                curvature_val =  1.0 / (u.armLength * tand(angle_val/ 2.0));
            end
           
            % Change sleeve curvature 
            positionServoSetting = getArmPositionFromAngle(angle_val);
            
            disp(["Angle: " num2str(angle_val) ", ..." + ...
                "Curvature:" num2str(curvature_val) ...
                ", motor:" num2str(positionServoSetting)] )

            setArmPosition(positionServoSetting);
            pause(2.0);
            
            recordForceData(true, i, 0, 50, 1, curvature_val, orientation);
        end
     end
    
    saveNewFile(dirName, num2str(orientation), 0, "_10"); 
    setArmPosition(1);
    pause(3)
end

function stopMotor()
    global m
    writePosition(m.force_servo_obj, 0.5);   
end

% moveForceToStartPosition - move load cell to start postion by 
%                            turning the lead screw until the load cell arm 
%                            hits to upper limit switch.
function moveForceToStartPosition()
    global m
    stopMotor();
    while(readDigitalPin(m.arduino_obj, m.upper_limit_pin))
        writePosition(m.force_servo_obj, 0.3);
    end
    stopMotor();
end

% resetVoltages - clear both graphs and the arrays containing voltages
function resetVoltages()
    global g d
    d.voltages = []
    
    set(g.force_graph,'XData',[0],'YData', [0]);
    
    for j = 1:5
        set(g.loss_graph(j), 'XData', [0], 'YData', [0]);
    end
    
    for j = 1:2
        set(g.predictedForce_graph(j),'XData',[0],'YData', [0]);
    end
    
    for j = 1:2
        set(g.predictedCurv_graph(j),'XData',[0],'YData', [0]);
    end

    for j = 1:2
        set(g.predictedOrien_graph(j),'XData',[0],'YData', [0]);
    end
end

% moveToContact - turn the lead screw until the load cell comes into
%                 contact with the sleeve
function moveToContact()
    global m 
    mass = 0;
    while(readDigitalPin(m.arduino_obj, m.lower_limit_pin) && mass < m.contact_mass)
        %mass = get_weight(m.calibration_obj, m.load_cell, 2);  
        mass = get_weight(m.calibration_obj, m.load_cell, 2);  
        writePosition(m.force_servo_obj, 0.7);
    end
    stopMotor();
end

% finalizePositions - if an error is encountered, use this function to 
%                     return to a neutral machine state. 
function finalizePositions()
    disp("Wait for cleanup")
    moveForceToStartPosition();
    setArmPosition(1.0);
    pause(3)
    setFiberOrientation(0);
end

% getLoadInNewtons - read the HX711 signal and convert to newtons
function load = getLoadInNewtons(nr_avg_readings)
    % load cell returns the mass in grams
    global m
    accel_of_gravity = 9.81;                                             % kg m / s^2
    mass = get_weight(m.calibration_obj, m.load_cell, nr_avg_readings);  % grams
    load = accel_of_gravity  * mass / 1000.0;                            % Newtons
end

% gotoLoad - turn the lead screw until the load reaches the maxLoad value
%            or the lower limit switch is hit. 
function gotoLoad(ID, increaseLoad, maxLoad, curvature, orientation)
    global m d
    
    reached_load      = false;     % flag 
    delay_in_sec      = 0.001;     % pause between data points in seconds
    max_motor_speed   = 0.08;      % number between 0 and 0.5 (larger is faster)

    if (increaseLoad)
        servo_setting = 0.5 + max_motor_speed;
        load_type     = 1;
    else
        servo_setting = 0.5 - max_motor_speed;
        load_type     = -1;
    end
    
    writePosition(m.force_servo_obj, servo_setting);

    while (~reached_load & readDigitalPin(m.arduino_obj, m.lower_limit_pin))
        
        % get and check mass reading
        load_cell_reading_N = recordForceData(true, ID, load_type, ...
            1, 2, curvature, orientation);

        if (increaseLoad)
            if (load_cell_reading_N > maxLoad)
                reached_load = true;
            end
        else 
            if (load_cell_reading_N < maxLoad)
                reached_load = true;
            end
        end

        pause(delay_in_sec);    
    end
    
    stopMotor();
end

% recordFroceData - take a reading from the load cell, update 
%                   graphs and store readings.
function force = recordForceData(plot_data, ID, load_type, ...
    nr_points_to_collect, nr_to_avg_per_reading, curvature, orientation)

    global m g d u

    for i = 1:nr_points_to_collect
        force = getLoadInNewtons(nr_to_avg_per_reading);  
        [v1, v2, v3, v4, v5]= readLossVoltages();
    
        % record data in array
        t = toc;
        lossArray  = - 10.0 * log10([v1, v2, v3, v4, v5] ./ d.baseline_values);

        % predicted force from neural network
        curvPredicted  = u.curv_neural_network.predictFcn([lossArray]);
        orienPredicted = u.orien_neural_network.predictFcn([lossArray]);
        
        if u.predictCurvOrien
            forcePredicted = ...
                u.force_neural_network.predictFcn([curvPredicted, orienPredicted, lossArray]);
        else 
            forcePredicted = ...
                u.force_neural_network.predictFcn([curvature, orientation, lossArray]);
        end

        % record data to data structure
        d.voltages = [d.voltages; v1, v2, v3, v4, v5, force, t, ...
            curvature, lossArray, ID, load_type, forcePredicted, ...
            curvPredicted, orienPredicted, orientation];
    
        % plot data
        if (plot_data)
            set(g.force_graph,'XData', d.voltages(:,7),'YData', d.voltages(:,6));
            
            for j = 1:5
                set(g.loss_graph(j), 'XData', d.voltages(:,7), ...
                    'YData', d.voltages(:, 8+j));
            end

            set(g.predictedForce_graph(1),'XData', d.voltages(:,7), ...
                'YData', d.voltages(:,6));
            set(g.predictedForce_graph(2),'XData', d.voltages(:,7), ...
                'YData', d.voltages(:, 16));

            
            set(g.predictedCurv_graph(1),'XData', d.voltages(:,7), ...
                'YData', d.voltages(:,8));
            set(g.predictedCurv_graph(2),'XData', d.voltages(:,7), ...
                'YData', d.voltages(:, 17));

            set(g.predictedOrien_graph(1),'XData', d.voltages(:,7), ...
                'YData', d.voltages(:, 19));
            set(g.predictedOrien_graph(2),'XData', d.voltages(:,7), ...
                'YData', d.voltages(:, 18));
        end
    end
end

% readLossVoltages - waveguides are connectived to the analog inputs on the
%                    arduino.
function [v1, v2, v3, v4, v5] = readLossVoltages()
    global m
    v1 = readVoltage(m.arduino_obj, 'A0');
    v2 = readVoltage(m.arduino_obj, 'A1');
    v3 = readVoltage(m.arduino_obj, 'A2');
    v4 = readVoltage(m.arduino_obj, 'A3');
    v5 = readVoltage(m.arduino_obj, 'A4');
end

function custom_error(message)
    stopMotor();
    error(message)
end

function saveNewFile(dirName, orientation, curvature_value,  save_type)
    global d 
    file_name = strcat("./", dirName, "/", orientation, ...
        '_curvature_', int2str(curvature_value), "_", save_type);
    
    % add losses
    output_data = d.voltages;
    for i = 1:5
        output_data(:,8+i) = - 10.0 * log10(output_data(:,i) / d.baseline_values(i));
    end

    % add titles to columns
    titles = ["A0[V]", "A1[V]", "A2[V]", "A3[V]", "A4[V]", "F[N]", "t[s]", ...
        "curvature[m^-1]", "loss0", "loss1","loss2","loss3","loss4", "ID", ...
        "loadType", "forcePredicted", "curvPredicted", "orienPredicted", ...
        "orientation"];
    loss_table = array2table(output_data, "VariableNames", titles);
    
    writetable(loss_table , strcat(file_name, ".xlsx"));

    % save copy of force diagram
    exportgraphics(gcf, strcat(file_name, ".png"), 'Resolution', 300);
    
    % save matlab data
    save(strcat(file_name, ".mat"), "d", "loss_table");

end

function saveBaselineData(dirName, orientation, curvature_value)
    global d
    file_name = strcat("./", dirName, "/", orientation, '_curvature_', ...
        int2str(curvature_value), "_", "baselines", ".xlsx");
    writematrix(d.baselineArray, file_name);
end

% getArmPositionFromAngle - convert the angleInDeg to a motor setting
function pos = getArmPositionFromAngle(angleInDeg)
    if (angleInDeg > 180)
        angleInDeg = 180;
    end
    if (angleInDeg < 90)
        angleInDeg = 90;
    end

    % calibration values determined from experiment
    m   = (0.98762 - 0.42573)/(180.0 - 90.0);
    b   = 0.42573;
    pos = m .* (angleInDeg - 90.0) + b;
end


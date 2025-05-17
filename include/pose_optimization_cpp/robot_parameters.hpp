#pragma once

#include <vector>

struct RobotParameters {
    double com2shoulder1_x = -0.027; 
    double com2shoulder1_y = 0.047; 
    double com2shoulder1_z = 0.162;
    double shoulder12shoulder2_x = 0.054; 
    double shoulder12shoulder2_y = 0.085; 
    double shoulder12shoulder2_z = 0;
    double shoulder22elbow_x = 0;  
    double shoulder22elbow_y = 0; 
    double shoulder22elbow_z = -0.25;
    double elbow2hand_x = 0; 
    double elbow2hand_y = 0; 
    double elbow2hand_z = -0.25; 
    
    double com2hip1_x = 0; 
    double com2hip1_y = 0.047; 
    double com2hip1_z = -0.1;
    double hip12hip2_x = 0.075; 
    double hip12hip2_y = 0.020; 
    double hip12hip2_z = -0.06;
    double hip22thigh_x = -0.09; 
    double hip22thigh_y = 0.023; 
    double hip22thigh_z = 0;
    double thigh2knee_x = 0; 
    double thigh2knee_y = 0; 
    double thigh2knee_z = -0.22;
    double knee2ankle_x = 0;  
    double knee2ankle_y = 0;  
    double knee2ankle_z = -0.22;
    double ankle2toe_x = 0.09; 
    double ankle2toe_y = 0; 
    double ankle2toe_z = -0.036;
    double ankle2heel_x = -0.05; 
    double ankle2heel_y = 0; 
    double ankle2heel_z = -0.036;
    double mass = 17.0;
};

struct ObjectParameters {
    double width = 0.3;
    double height = 0.5;
    double depth = 0.3;
};



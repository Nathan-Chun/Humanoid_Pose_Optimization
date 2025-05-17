// optimizer.cpp
#include "optimizer.hpp"
#include "kinematics.hpp"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <casadi/casadi.hpp>


using namespace casadi;

Optimizer::Optimizer() {
    setup_problem();
}

void Optimizer::setup_problem() {
    // Define variables

    // Define constants
    double g = 9.81; // Gravity
    double mu_foot = 1.0;  // Foot-ground friction
    double mu_hand = 1.0;  // Hand-object friction
    double mu_object = 0.5; // Object-ground friction

    // p_c - COM Position
    p_c = opti.variable(3, 1);
    MX p_c_x = p_c(0); MX p_c_y = p_c(1); MX p_c_z = p_c(2);

    // rpy - COM Orientation
    rpy = opti.variable(3, 1);
    MX roll = rpy(0); MX pitch = rpy(1); MX yaw = rpy(2);

    // arm joint angles
    // q_arms = opti.variable(8, 1);
    q_arms = opti.variable(6, 1);
    MX q1L_a = q_arms(0); MX q2L_a = q_arms(1); MX q3L_a = q_arms(2); //MX q4L_a = q_arms(3); 
    MX q1R_a = q_arms(3); MX q2R_a = q_arms(4); MX q3R_a = q_arms(5); //MX q4R_a = q_arms(7); 

    // leg joint angles
    q_legs = opti.variable(10, 1);
    MX q1L_l = q_legs(0); MX q2L_l = q_legs(1); MX q3L_l = q_legs(2); MX q4L_l = q_legs(3); MX q5L_l = q_legs(4);
    MX q1R_l = q_legs(5); MX q2R_l = q_legs(6); MX q3R_l = q_legs(7); MX q4R_l = q_legs(8); MX q5R_l = q_legs(9); 

    // F_arms = [Falx; Faly; Falz; Farx; Fary, Farz]
    F_arms = opti.variable(6, 1);
    MX Falx = F_arms(0); MX Faly = F_arms(1); MX Falz = F_arms(2); 
    MX Farx = F_arms(3); MX Fary = F_arms(4); MX Farz = F_arms(5);
    // MX Farx = F_arms(0); MX Fary = F_arms(1); MX Farz = F_arms(2);

    // F_legs = [Fllx; Flly; Fllz; Flrx; Flry, Flrz]
    F_legs = opti.variable(6, 1);
    MX Fllx = F_legs(0); MX Flly = F_legs(1); MX Fllz = F_legs(2);
    MX Flrx = F_legs(3); MX Flry = F_legs(4); MX Flrz = F_legs(5);
    // MX Flrx = F_legs(0); MX Flry = F_legs(1); MX Flrz = F_legs(2);
    //-----------------------------------------------------------------------------------------------------//
    //Problem Definition
    // Desired COM Positon
    MX p_des = MX::vertcat({0.0, 0.0, 0.5}); // Example desired position

    // Object info
    struct Object {
        Eigen::Vector3d size;
        Eigen::Vector3d com;
        double mass;          
        Eigen::Matrix3d I;    // Inertia matrix
    };

    // Define robot properties
    Object object;
    object.size = Eigen::Vector3d(1.0, 1.0, 1.0);
    object.com = Eigen::Vector3d(0.8, 0.0, object.size(2) / 2.0);
    object.mass = 20.0;

    // Calculate inertia matrix
    object.I = Eigen::Matrix3d::Zero();
    object.I(0, 0) = (1.0 / 12.0) * object.mass * (object.size(1) * object.size(1) + object.size(2) * object.size(2));
    object.I(1, 1) = (1.0 / 12.0) * object.mass * (object.size(0) * object.size(0) + object.size(2) * object.size(2));
    object.I(2, 2) = (1.0 / 12.0) * object.mass * (object.size(0) * object.size(0) + object.size(1) * object.size(1));

    // Print results
    std::cout << "Robot COM: " << object.com.transpose() << std::endl;
    std::cout << "Robot Inertia Matrix:\n" << object.I << std::endl;

    //-----------------------------------------------------------------------------------------------------//
    //Constraints
    // MX locations = Kinematics::get_locations({p_c_x; p_c_y; p_c_z}, {roll; pitch; yaw}, {q1L_a;q2L_a;q3L_a;q1R_a;q2R_a;q3R_a}, {q1L_l;q2L_l;q3L_l;q4L_l;q5L_l; q1R_l;q2R_l;q3R_l;q4R_l;q5R_l});
    RobotParameters params;
    auto locations = Kinematics::get_locations(p_c, rpy, q_arms, q_legs, params);
    
    std::cout << "Locations: " << locations.size()<< std::endl;

    // foot location on the ground:
    opti.subject_to(locations[29]==0); opti.subject_to(locations[32]==0); 
    opti.subject_to(locations[62]==0); opti.subject_to(locations[65]==0);
    // MX footL = (locations(28:30) +locations(31:33))/2;
    // MX footR = (locations(61:63) +locations(64:66))/2;

    // Foot location on the ground
    MX footL = (MX::vertcat({locations[27], locations[28], locations[29]}) +
            MX::vertcat({locations[30], locations[31], locations[32]})) / 2;

    MX footR = (MX::vertcat({locations[60], locations[61], locations[62]}) +
            MX::vertcat({locations[63], locations[64], locations[65]})) / 2;

    // std::cout << "FootL: " << footL << std::endl;

    //CoM at desired task CoM location
    opti.subject_to(p_c_y == p_des(1));

    //Hand on the object:
    MX handL = MX::vertcat({locations[9], locations[10], locations[11]});
    MX handR = MX::vertcat({locations[42], locations[43], locations[44]});

    // Case 1: Hands on the vertical side:
    // MX il = if_else(handL(1) <=(object.com(1) - object.size(1)/2), 1, 0);
    MX i1 = 1; MX i3 = 1;
    opti.subject_to(handL(0)*i1 == (object.com(0) - object.size(0)/2)*i1);
    opti.subject_to(handL(2)*i1 <= object.size(2)*i1);

    // MX i3 == if_else(handR(1) <= object.com(1) - object.size(1)/2, 1, 0);
    opti.subject_to(handR(0)*i3 == (object.com(0) - object.size(0)/2)*i3);
    opti.subject_to(handR(2)*i3 <= object.size(2)*i3);
    opti.subject_to( -object.size(1)/2.3 + object.com(1) <= handR(1) <= object.size(1)/2.3 + object.com(1));
    opti.subject_to( -object.size(1)/2.3 + object.com(1) <= handL(1) <= object.size(1)/2.3 + object.com(1));
    
    // Statics constraints
    // force equilibrium
    opti.subject_to( Fllx + Flrx + Falx + Farx == 0);
    opti.subject_to( Flly + Flry + Faly + Fary == 0);
    opti.subject_to( Fllz + Flrz + Falz + Farz - 17*g == 0);

    std::cout << "Static constraints set" << std::endl;

    // Moment equilibrium
    opti.subject_to(
        skew(handL - MX::vertcat({p_c_x, p_c_y, p_c_z})) * MX::vertcat({Falx, Faly, Falz}) +
        skew(handR - MX::vertcat({p_c_x, p_c_y, p_c_z})) * MX::vertcat({Farx, Fary, Farz}) +
        skew(footL - MX::vertcat({p_c_x, p_c_y, p_c_z})) * MX::vertcat({Fllx, Flly, Fllz}) +
        skew(footR - MX::vertcat({p_c_x, p_c_y, p_c_z})) * MX::vertcat({Flrx, Flry, Flrz}) ==
        MX::vertcat({0, 0, 0})
    );
    std::cout << "Moment constraints set" << std::endl;

    // Friction constraints
    opti.subject_to( abs(Flrx) <= mu_foot*abs(Flrz) );
    opti.subject_to( abs(Flry) <= mu_foot*abs(Flrz) );
    opti.subject_to( abs(Fllx) <= mu_foot*abs(Fllz) );
    opti.subject_to( abs(Flly) <= mu_foot*abs(Fllz) );
    opti.subject_to( abs(Farz) <= mu_hand*abs(Farx) );
    opti.subject_to( abs(Fary) <= mu_hand*abs(Farx) );
    opti.subject_to( abs(Falz) <= mu_hand*abs(Falx) );
    opti.subject_to( abs(Faly) <= mu_hand*abs(Falx) );

    opti.subject_to( Falx == -mu_object * object.mass/2 * g );
    opti.subject_to( Farx == -mu_object * object.mass/2 * g );

    std::cout << "Friction constraints set" << std::endl;
    
    MX p_hand_object_l = MX::vertcat({object.com(0), object.com(1), object.com(2)}) - handL;
    MX p_hand_object_r = MX::vertcat({object.com(0), object.com(1), object.com(2)}) - handR;
    // opti.subject_to( Faly == Falx*(p_hand_object_l(2)/p_hand_object_l(1)) );
    // opti.subject_to( Fary == Farx*(p_hand_object_r(2)/p_hand_object_r(1)) );
    opti.subject_to( Falz == Falx*(p_hand_object_l(2)/p_hand_object_l(0)) );
    opti.subject_to( Farz == Farx*(p_hand_object_r(2)/p_hand_object_r(0)) );

    // foot at origin:
    opti.subject_to( footL(0) == 0);
    opti.subject_to( footR(0) == 0);

    // shoulder close to object
    MX shoulderL = MX::vertcat({locations[3], locations[4], locations[5]});
    MX shoulderR = MX::vertcat({locations[36], locations[37], locations[38]});

    opti.subject_to( 0.07 <= handL(0)-shoulderL(0) <= 0.15 );
    opti.subject_to( 0.07 <= handR(0)-shoulderR(0) <= 0.15 );

    // hip close to foot
    MX thighL = MX::vertcat({locations[18], locations[19], locations[20]});
    MX thighR = MX::vertcat({locations[51], locations[52], locations[53]});
    
    opti.subject_to( 0.0 <= thighL(0)-footL(0) <= 0.04 );
    opti.subject_to( 0.0 <= thighR(0)-footR(0) <= 0.04 );

    std::cout << "p_c: " << p_c << std::endl;

    // Helper lambda function for degrees to radians
    auto deg2rad = [](double deg) { return deg * M_PI / 180.0; };

    // Arm joint limits
    opti.subject_to( deg2rad(10)   <= q1L_a <= deg2rad(90) );
    opti.subject_to( deg2rad(-90)  <= q1R_a <= deg2rad(-10) );
    opti.subject_to( deg2rad(-90)  <= q2L_a <= deg2rad(90) );
    opti.subject_to( deg2rad(-90)  <= q2R_a <= deg2rad(90) );
    opti.subject_to( deg2rad(-160) <= q3L_a <= deg2rad(-10) );
    opti.subject_to( deg2rad(-160) <= q3R_a <= deg2rad(-10) );

    // Leg joint limits
    opti.subject_to( deg2rad(-30)  <= q1L_l <= deg2rad(30) );
    opti.subject_to( deg2rad(-30)  <= q1R_l <= deg2rad(30) );
    opti.subject_to( deg2rad(-30)  <= q2L_l <= deg2rad(30) );
    opti.subject_to( deg2rad(-30)  <= q2R_l <= deg2rad(30) );
    opti.subject_to( deg2rad(-10)  <= q3L_l <= deg2rad(90) );
    opti.subject_to( deg2rad(-10)  <= q3R_l <= deg2rad(90) );
    opti.subject_to( deg2rad(-150) <= q4L_l <= deg2rad(-10) );
    opti.subject_to( deg2rad(-150) <= q4R_l <= deg2rad(-10) );
    opti.subject_to( deg2rad(-60)  <= q5L_l <= deg2rad(60) );
    // opti.subject_to( deg2rad(-60)  <= q5R_l <= deg2rad(60) );


    //-----------------------------------------------------------------------------------------------------//
    // Solver setup
    casadi::Dict p_opts;
    p_opts["expand"] = true;
    // p_opts.discrete = [false,false,false,false,false,true];
    opti.solver("ipopt", p_opts);

    //-----------------------------------------------------------------------------------------------------//
    // Initial Guess
    DM p_guess = DM::vertcat({0.0, 0.0, 0.5});
    DM rpy_guess = DM::vertcat({0,pi/4,0});
    DM q_arms_guess = DM::vertcat({0, 0,-pi/2, 0, 0, -pi/2});
    DM q_legs_guess = DM::vertcat({0,0,pi/4,-pi/2,pi/4, 0,0,pi/4,-pi/2,pi/4});
    DM F_arms_guess = DM::vertcat({-mu_object * object.mass/2 * g, 0, 0, -mu_object * object.mass/2 * g,0,0});
    //New
    DM F_legs_guess = DM::zeros(6, 1);
    opti.set_initial(F_legs, F_legs_guess);

    opti.set_initial(F_arms,F_arms_guess);
    opti.set_initial(p_c,p_guess);
    opti.set_initial(rpy,rpy_guess);
    opti.set_initial(q_arms,q_arms_guess);
    opti.set_initial(q_legs,q_legs_guess);

    std::cout << "Initial guess: " << std::endl;
    //-----------------------------------------------------------------------------------------------------//
 
    // Cost Function
    MX J_pitch = pow(pitch,2) + pow(roll,2) + pow(yaw,2); // minimize rotation angle
    MX J_foot = pow(p_c_x*0 - (locations[27]+locations[30])/2, 2) + pow(p_c_x*0 - (locations[60]+locations[63])/2, 2); // ensure stability region
    MX J_ig = sum1(pow(q_arms - q_arms_guess, 2)) + sum1(pow(q_legs - q_legs_guess, 2));
    // MX J_ig = pow((q_arms - q_arms_guess).T() * (q_arms - q_arms_guess) + (q_legs - q_legs_guess).T() * (q_legs - q_legs_guess), 1); // already quadratic
    MX J_com_z = pow(p_c_z-p_des(2), 2);
    MX J_F_l = pow(Fllx,2) + pow(Flrx,2) + pow(Flly,2) + pow(Flry,2) + pow(Fllz,2) + pow(Flrz,2);
    MX J_F_a = pow(Falx,2) + pow(Farx,2) + pow(Faly,2) + pow(Fary,2) + pow(Falz,2) + pow(Farz,2);
    MX J_hand_pos = pow(handL(2)-object.com(2),2) + pow(handR(2)-object.com(2),2);
    MX J_CoM_pos = pow((object.com(0) - p_c_x)*10*object.mass, 2);
    
    MX J = 1*J_pitch + 0*J_foot + 100*J_ig + 100*J_com_z + 0*J_F_l + 10*J_F_a + 1*J_hand_pos + 0*J_CoM_pos;
    // J = 0;
    // std::cout << "Objective function: " << J << std::endl;
    // std::cout << "Initial p_c: " << opti.debug().value(p_c) << std::endl;
    // std::cout << "Initial q_arms: " << opti.debug().value(q_arms) << std::endl;
    // std::cout << "Initial q_legs: " << opti.debug().value(q_legs) << std::endl;
    // std::cout << "Initial handL: " << opti.debug().value(handL) << std::endl;
    // std::cout << "Initial footL: " << opti.debug().value(footL) << std::endl;
    opti.minimize(J);
    
    //-----------------------------------------------------------------------------------------------------//
    // Solution
    // OptiSol sol = opti.solve();

    // DM p_c_val = sol.value(p_c);
    // std::vector<double> p_c_sol = {
    //     p_c_val(0).scalar(), p_c_val(1).scalar(), p_c_val(2).scalar()
    // };
    
    // DM rpy_val = sol.value(rpy);
    // std::vector<double> rpy_sol = {
    //     rpy_val(0).scalar(), rpy_val(1).scalar(), rpy_val(2).scalar()
    // };
    
    // DM q_arms_val = sol.value(q_arms);
    // std::vector<double> q_arms_sol = {
    //     q_arms_val(0).scalar(), q_arms_val(1).scalar(), q_arms_val(2).scalar(), q_arms_val(3).scalar(),
    //     q_arms_val(4).scalar(), q_arms_val(5).scalar(), q_arms_val(6).scalar(), q_arms_val(7).scalar()
    // };

    // DM q_legs_val = sol.value(q_legs);
    // std::vector<double> q_legs_sol = {
    //     q_legs_val(0).scalar(), q_legs_val(1).scalar(), q_legs_val(2).scalar(), q_legs_val(3).scalar(),
    //     q_legs_val(4).scalar(), q_legs_val(5).scalar(), q_legs_val(6).scalar(), q_legs_val(7).scalar(),
    //     q_legs_val(8).scalar(), q_legs_val(9).scalar()
    // };
    
    // DM F_arms_val = sol.value(F_arms);
    // std::vector<double> F_arms_sol = {
    //     F_arms_val(0).scalar(), F_arms_val(1).scalar(), F_arms_val(2).scalar(),
    //     F_arms_val(3).scalar(), F_arms_val(4).scalar(), F_arms_val(5).scalar()
    // };
    
    // DM F_legs_val = sol.value(F_legs);
    // std::vector<double> F_legs_sol = {
    //     F_legs_val(0).scalar(), F_legs_val(1).scalar(), F_legs_val(2).scalar(),
    //     F_legs_val(3).scalar(), F_legs_val(4).scalar(), F_legs_val(5).scalar()
    // };

    //-----------------------------------------------------------------------------------------------------//

    // Visual_Humanoid_Pose_Optimization(p_c_sol,rpy_sol,q_arms_sol,q_legs_sol);

    // // Add constraints
    // opti.subject_to(p_torso(2) >= 0.5);

    // // Objective
    // objective = pow(q(0), 2); // dummy
    // opti.minimize(objective);

    // // Set solver
    // opti.solver("ipopt");
}

casadi::MX Optimizer::skew(const casadi::MX& a) {
    return MX::vertcat({
        MX::horzcat({0, -a(2), a(1)}),
        MX::horzcat({a(2), 0, -a(0)}),
        MX::horzcat({-a(1), a(0), 0})
    });
}

void Optimizer::solve() {

    auto solution = opti.solve();
    // auto stats = solution.get_stats();
    // std::cout << "Solution status: " << stats.at("return_status") << std::endl;
    std::cout << "Torso position: " << solution.value(p_c) << std::endl;
    std::cout << "Torso orientation: " << solution.value(rpy) << std::endl;
    std::cout << "Arm joint angles: " << solution.value(q_arms) << std::endl;
    std::cout << "Leg joint angles: " << solution.value(q_legs) << std::endl;
    std::cout << "Arm forces: " << solution.value(F_arms) << std::endl;
    std::cout << "Leg forces: " << solution.value(F_legs) << std::endl;
}

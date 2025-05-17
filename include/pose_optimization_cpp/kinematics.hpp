// kinematics.hpp
#pragma once

#include <casadi/casadi.hpp>
#include <vector>
#include "robot_parameters.hpp"

using casadi::MX;

class Kinematics {
    public:
        static std::vector<casadi::MX> get_locations(
            const casadi::MX& p_c, 
            const casadi::MX& rpy, 
            const casadi::MX& q_arms, 
            const casadi::MX& q_legs,
            const RobotParameters& params
        );
    };
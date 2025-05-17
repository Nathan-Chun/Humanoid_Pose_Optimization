#pragma once

# include <casadi/casadi.hpp>
# include "robot_parameters.hpp"

class Optimizer{
    public:
        Optimizer();
        void setup_problem();
        void solve();
        casadi::MX skew(const casadi::MX& a); // Declare skew as a member function


    private:
        casadi::Opti opti;
        casadi::MX p_c, rpy, q_arms, q_legs, F_arms, F_legs;
        RobotParameters robot_params;
        casadi::MX decision_vars;
        casadi:: MX objective;
};
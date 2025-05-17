// kinematics.cpp
#include "kinematics.hpp"
#include "robot_parameters.hpp"
using namespace casadi;

std::vector<casadi::MX> Kinematics::get_locations(const casadi::MX& p_c, const casadi::MX& rpy, const casadi::MX& q_arms, const casadi::MX& q_legs, const RobotParameters& params) {
    // RobotParameters params;
    // Compute forward kinematics here
    // Return a vector of MX expressions representing link positions

    // Initialize locations vector
    std::vector<MX> locations;

    MX roll = rpy(0); MX pitch = rpy(1); MX yaw = rpy(2);
    MX q1al = q_arms(0); MX q2al = q_arms(1); MX q3al = q_arms(2); //MX q4al = q_arms(3);
    MX q1ar = q_arms(3); MX q2ar = q_arms(4); MX q3ar = q_arms(5); //MX q4ar = q_arms(7);

    MX q1ll = q_legs(0); MX q2ll = q_legs(1); MX q3ll = q_legs(2); MX q4ll = q_legs(3); MX q5ll = q_legs(4);
    MX q1lr = q_legs(5); MX q2lr = q_legs(6); MX q3lr = q_legs(7); MX q4lr = q_legs(8); MX q5lr = q_legs(9);

    MX x = p_c(0); MX y = p_c(1); MX z = p_c(2);

    int side = 1;

    // Rotation matrices Rx, Ry, Rz
    auto Rx = [](const MX& a) {
        return MX::vertcat({
            MX::horzcat({1, 0, 0}),
            MX::horzcat({0, cos(a), -sin(a)}),
            MX::horzcat({0, sin(a), cos(a)})
        });
    };
    auto Ry = [](const MX& a) {
        return MX::vertcat({
            MX::horzcat({cos(a), 0, sin(a)}),
            MX::horzcat({0, 1, 0}),
            MX::horzcat({-sin(a), 0, cos(a)})
        });
    };
    auto Rz = [](const MX& a) {
        return MX::vertcat({
            MX::horzcat({cos(a), -sin(a), 0}),
            MX::horzcat({sin(a), cos(a), 0}),
            MX::horzcat({0, 0, 1})
        });
    };

    // Homogeneous transformation helpers
    auto homog = [](const MX& R, const MX& t) {
        return MX::vertcat({
            MX::horzcat({R, t}),
            MX::horzcat({MX::zeros(1,3), MX::ones(1,1)})
        });
    };

    MX H_com = homog(MX::eye(3), MX::vertcat({x, y, z}));

    // handL
    MX H_com2shoulder1 = MX::mtimes({homog(MX::mtimes({Rz(yaw), Ry(pitch),Rx(roll)}), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.com2shoulder1_x, params.com2shoulder1_y*side, params.com2shoulder1_z}))});

    MX H_shoulder12shoulder2 = MX::mtimes({
        homog(Rx(q1al), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.shoulder12shoulder2_x, params.shoulder12shoulder2_y*side, params.shoulder12shoulder2_z}))
    });

    MX H_shoulder22elbow = MX::mtimes({
        homog(Ry(q2al), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.shoulder22elbow_x, params.shoulder22elbow_y*side, params.shoulder22elbow_z}))
    });

    MX H_elbow2hand = MX::mtimes({
        homog(Ry(q3al), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.elbow2hand_x, params.elbow2hand_y*side, params.elbow2hand_z}))
    });

    // footL
    MX H_com2hip1 = MX::mtimes({
        homog(MX::mtimes({Rz(yaw), Ry(pitch), Rx(roll)}), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.com2hip1_x, params.com2hip1_y*side, params.com2hip1_z}))
    });

    MX H_hip12hip2 = MX::mtimes({
        homog(Rz(q1ll), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.hip12hip2_x, params.hip12hip2_y*side, params.hip12hip2_z}))
    });

    MX H_hip22thigh = MX::mtimes({
        homog(Rx(q2ll), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.hip22thigh_x, params.hip22thigh_y*side, params.hip22thigh_z}))
    });

    MX H_thigh2knee = MX::mtimes({
        homog(Ry(q3ll), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.thigh2knee_x, params.thigh2knee_y*side, params.thigh2knee_z}))
    });

    MX H_knee2ankle = MX::mtimes({
        homog(Ry(q4ll), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.knee2ankle_x, params.knee2ankle_y*side, params.knee2ankle_z}))
    });

    MX H_ankle2toe = MX::mtimes({
        homog(Ry(q5ll), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.ankle2toe_x, params.ankle2toe_y, params.ankle2toe_z}))
    });

    MX H_ankle2heel = MX::mtimes({
        homog(Ry(q5ll), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.ankle2heel_x, params.ankle2heel_y, params.ankle2heel_z}))
    });

    MX point = MX::vertcat({0, 0, 0, 1});

    MX shoulder1l = MX::mtimes({H_com, H_com2shoulder1, point});
    MX shoulder2l = MX::mtimes({H_com, H_com2shoulder1, H_shoulder12shoulder2, point});
    MX elbowl     = MX::mtimes({H_com, H_com2shoulder1, H_shoulder12shoulder2, H_shoulder22elbow, point});
    MX handl      = MX::mtimes({H_com, H_com2shoulder1, H_shoulder12shoulder2, H_shoulder22elbow, H_elbow2hand, point});
    MX hip1l      = MX::mtimes({H_com, H_com2hip1, point});
    MX hip2l      = MX::mtimes({H_com, H_com2hip1, H_hip12hip2, point});
    MX thighl     = MX::mtimes({H_com, H_com2hip1, H_hip12hip2, H_hip22thigh, point});
    MX kneel      = MX::mtimes({H_com, H_com2hip1, H_hip12hip2, H_hip22thigh, H_thigh2knee, point});
    MX anklel     = MX::mtimes({H_com, H_com2hip1, H_hip12hip2, H_hip22thigh, H_thigh2knee, H_knee2ankle, point});
    MX toel       = MX::mtimes({H_com, H_com2hip1, H_hip12hip2, H_hip22thigh, H_thigh2knee, H_knee2ankle, H_ankle2toe, point});
    MX heell      = MX::mtimes({H_com, H_com2hip1, H_hip12hip2, H_hip22thigh, H_thigh2knee, H_knee2ankle, H_ankle2heel, point});

    side = -1;

    // handR
    MX H_com2shoulder1_r = MX::mtimes({
        homog(MX::mtimes({Rz(yaw), Ry(pitch), Rx(roll)}), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.com2shoulder1_x, params.com2shoulder1_y*side, params.com2shoulder1_z}))
    });

    MX H_shoulder12shoulder2_r = MX::mtimes({
        homog(Rx(q1ar), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.shoulder12shoulder2_x, params.shoulder12shoulder2_y*side, params.shoulder12shoulder2_z}))
    });

    MX H_shoulder22elbow_r = MX::mtimes({
        homog(Ry(q2ar), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.shoulder22elbow_x, params.shoulder22elbow_y*side, params.shoulder22elbow_z}))
    });

    MX H_elbow2hand_r = MX::mtimes({
        homog(Ry(q3ar), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.elbow2hand_x, params.elbow2hand_y*side, params.elbow2hand_z}))
    });

    // footR
    MX H_com2hip1_r = MX::mtimes({
        homog(MX::mtimes({Rz(yaw), Ry(pitch), Rx(roll)}), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.com2hip1_x, params.com2hip1_y*side, params.com2hip1_z}))
    });

    MX H_hip12hip2_r = MX::mtimes({
        homog(Rz(q1lr), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.hip12hip2_x, params.hip12hip2_y*side, params.hip12hip2_z}))
    });

    MX H_hip22thigh_r = MX::mtimes({
        homog(Rx(q2lr), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.hip22thigh_x, params.hip22thigh_y*side, params.hip22thigh_z}))
    });

    MX H_thigh2knee_r = MX::mtimes({
        homog(Ry(q3lr), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.thigh2knee_x, params.thigh2knee_y*side, params.thigh2knee_z}))
    });

    MX H_knee2ankle_r = MX::mtimes({
        homog(Ry(q4lr), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.knee2ankle_x, params.knee2ankle_y*side, params.knee2ankle_z}))
    });

    MX H_ankle2toe_r = MX::mtimes({
        homog(Ry(q5lr), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.ankle2toe_x, params.ankle2toe_y, params.ankle2toe_z}))
    });

    MX H_ankle2heel_r = MX::mtimes({
        homog(Ry(q5lr), MX::vertcat({0, 0, 0})),
        homog(MX::eye(3), MX::vertcat({params.ankle2heel_x, params.ankle2heel_y, params.ankle2heel_z}))
    });

    MX shoulder1r = MX::mtimes({H_com, H_com2shoulder1_r, point});
    MX shoulder2r = MX::mtimes({H_com, H_com2shoulder1_r, H_shoulder12shoulder2_r, point});
    MX elbowr     = MX::mtimes({H_com, H_com2shoulder1_r, H_shoulder12shoulder2_r, H_shoulder22elbow_r, point});
    MX handr      = MX::mtimes({H_com, H_com2shoulder1_r, H_shoulder12shoulder2_r, H_shoulder22elbow_r, H_elbow2hand_r, point});
    MX hip1r      = MX::mtimes({H_com, H_com2hip1_r, point});
    MX hip2r      = MX::mtimes({H_com, H_com2hip1_r, H_hip12hip2_r, point});
    MX thighr     = MX::mtimes({H_com, H_com2hip1_r, H_hip12hip2_r, H_hip22thigh_r, point});
    MX kneer      = MX::mtimes({H_com, H_com2hip1_r, H_hip12hip2_r, H_hip22thigh_r, H_thigh2knee_r, point});
    MX ankler     = MX::mtimes({H_com, H_com2hip1_r, H_hip12hip2_r, H_hip22thigh_r, H_thigh2knee_r, H_knee2ankle_r, point});
    MX toer       = MX::mtimes({H_com, H_com2hip1_r, H_hip12hip2_r, H_hip22thigh_r, H_thigh2knee_r, H_knee2ankle_r, H_ankle2toe_r, point});
    MX heelr      = MX::mtimes({H_com, H_com2hip1_r, H_hip12hip2_r, H_hip22thigh_r, H_thigh2knee_r, H_knee2ankle_r, H_ankle2heel_r, point});


    // Add all locations to locations vector
    locations.push_back(shoulder1l(Slice(0,3)));
    locations.push_back(shoulder2l(Slice(0,3)));
    locations.push_back(elbowl(Slice(0,3)));
    locations.push_back(handl(Slice(0,3)));
    locations.push_back(hip1l(Slice(0,3)));
    locations.push_back(hip2l(Slice(0,3)));
    locations.push_back(thighl(Slice(0,3)));
    locations.push_back(kneel(Slice(0,3)));
    locations.push_back(anklel(Slice(0,3)));
    locations.push_back(toel(Slice(0,3)));
    locations.push_back(heell(Slice(0,3)));

    locations.push_back(shoulder1r(Slice(0,3)));
    locations.push_back(shoulder2r(Slice(0,3)));
    locations.push_back(elbowr(Slice(0,3)));
    locations.push_back(handr(Slice(0,3)));
    locations.push_back(hip1r(Slice(0,3)));
    locations.push_back(hip2r(Slice(0,3)));
    locations.push_back(thighr(Slice(0,3)));
    locations.push_back(kneer(Slice(0,3)));
    locations.push_back(ankler(Slice(0,3)));
    locations.push_back(toer(Slice(0,3)));
    locations.push_back(heelr(Slice(0,3)));

    std::vector<casadi::MX> flat_locations;
    for (auto& loc : locations) {
        for (int i = 0; i < 3; ++i) {
            flat_locations.push_back(loc(i));
        }
    }

    // Example: locations.push_back(torso_pose + joint_angles(Slice(0, 3)));
    return flat_locations;
}

// Copyright 2022-2024 The MathWorks, Inc.
// Generated 27-Oct-2024 09:47:07
#ifndef _SLROS2_INITIALIZE_H_
#define _SLROS2_INITIALIZE_H_
#include "nmpcRect_types.h"
// Generic pub-sub header
#include "slros2_generic_pubsub.h"
#include "slros2_generic_param.h"
#ifndef SET_QOS_VALUES
#define SET_QOS_VALUES(qosStruct, _history, _depth, _durability, _reliability, _deadline \
, _lifespan, _liveliness, _lease_duration, _avoid_ros_namespace_conventions)             \
    {                                                                                    \
        qosStruct.history = _history;                                                    \
        qosStruct.depth = _depth;                                                        \
        qosStruct.durability = _durability;                                              \
        qosStruct.reliability = _reliability;                                            \
        qosStruct.deadline.sec = _deadline.sec;                                          \
        qosStruct.deadline.nsec = _deadline.nsec;                                        \
        qosStruct.lifespan.sec = _lifespan.sec;                                          \
        qosStruct.lifespan.nsec = _lifespan.nsec;                                        \
        qosStruct.liveliness = _liveliness;                                              \
        qosStruct.liveliness_lease_duration.sec = _lease_duration.sec;                   \
        qosStruct.liveliness_lease_duration.nsec = _lease_duration.nsec;                 \
        qosStruct.avoid_ros_namespace_conventions = _avoid_ros_namespace_conventions;    \
    }
#endif
inline rclcpp::QoS getQOSSettingsFromRMW(const rmw_qos_profile_t& qosProfile) {
    rclcpp::QoS qos(rclcpp::QoSInitialization::from_rmw(qosProfile));
    if (RMW_QOS_POLICY_DURABILITY_TRANSIENT_LOCAL == qosProfile.durability) {
        qos.transient_local();
    } else {
        qos.durability_volatile();
    }
    if (RMW_QOS_POLICY_RELIABILITY_RELIABLE == qosProfile.reliability) {
        qos.reliable();
    } else {
        qos.best_effort();
    }
    return qos;
}
// nmpcRect/Command Velocity Publisher/Publish
extern SimulinkPublisher<geometry_msgs::msg::Twist,SL_Bus_geometry_msgs_Twist> Pub_nmpcRect_168;
// nmpcRect/NMPC Controller/Optimization Status Publisher/Publish
extern SimulinkPublisher<std_msgs::msg::Int32,SL_Bus_std_msgs_Int32> Pub_nmpcRect_264;
// nmpcRect/NMPC Controller/Ref path publisher/Publish
extern SimulinkPublisher<geometry_msgs::msg::Point,SL_Bus_geometry_msgs_Point> Pub_nmpcRect_271;
// nmpcRect/NMPC Controller/State Publisher/Publish
extern SimulinkPublisher<catarob_interfaces::msg::StateEstimate,SL_Bus_catarob_interfaces_StateEstimate> Pub_nmpcRect_253;
// nmpcRect/Subscribe
extern SimulinkSubscriber<sensor_msgs::msg::NavSatFix,SL_Bus_sensor_msgs_NavSatFix> Sub_nmpcRect_272;
// nmpcRect/Subscribe1
extern SimulinkSubscriber<std_msgs::msg::Float64,SL_Bus_std_msgs_Float64> Sub_nmpcRect_273;
// nmpcRect/Subscribe2
extern SimulinkSubscriber<sensor_msgs::msg::Imu,SL_Bus_sensor_msgs_Imu> Sub_nmpcRect_274;
// For Block nmpcRect/Get Parameter
extern SimulinkParameterArrayGetter<real64_T,std::vector<double>> ParamGet_nmpcRect_289;
#endif

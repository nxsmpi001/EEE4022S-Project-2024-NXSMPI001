// Copyright 2022-2024 The MathWorks, Inc.
// Generated 27-Oct-2024 09:42:58
#include "slros2_initialize.h"
// nmpcSine/Command Velocity Publisher/Publish
SimulinkPublisher<geometry_msgs::msg::Twist,SL_Bus_geometry_msgs_Twist> Pub_nmpcSine_168;
// nmpcSine/NMPC Controller/Optimization Status Publisher/Publish
SimulinkPublisher<std_msgs::msg::Int32,SL_Bus_std_msgs_Int32> Pub_nmpcSine_264;
// nmpcSine/NMPC Controller/Ref path publisher/Publish
SimulinkPublisher<geometry_msgs::msg::Point,SL_Bus_geometry_msgs_Point> Pub_nmpcSine_305;
// nmpcSine/NMPC Controller/State Publisher/Publish
SimulinkPublisher<catarob_interfaces::msg::StateEstimate,SL_Bus_catarob_interfaces_StateEstimate> Pub_nmpcSine_253;
// nmpcSine/Subscribe
SimulinkSubscriber<sensor_msgs::msg::NavSatFix,SL_Bus_sensor_msgs_NavSatFix> Sub_nmpcSine_272;
// nmpcSine/Subscribe1
SimulinkSubscriber<std_msgs::msg::Float64,SL_Bus_std_msgs_Float64> Sub_nmpcSine_273;
// nmpcSine/Subscribe2
SimulinkSubscriber<sensor_msgs::msg::Imu,SL_Bus_sensor_msgs_Imu> Sub_nmpcSine_274;
// For Block nmpcSine/Get Parameter
SimulinkParameterArrayGetter<real64_T,std::vector<double>> ParamGet_nmpcSine_289;

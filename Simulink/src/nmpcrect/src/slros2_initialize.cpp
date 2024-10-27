// Copyright 2022-2024 The MathWorks, Inc.
// Generated 27-Oct-2024 09:47:07
#include "slros2_initialize.h"
// nmpcRect/Command Velocity Publisher/Publish
SimulinkPublisher<geometry_msgs::msg::Twist,SL_Bus_geometry_msgs_Twist> Pub_nmpcRect_168;
// nmpcRect/NMPC Controller/Optimization Status Publisher/Publish
SimulinkPublisher<std_msgs::msg::Int32,SL_Bus_std_msgs_Int32> Pub_nmpcRect_264;
// nmpcRect/NMPC Controller/Ref path publisher/Publish
SimulinkPublisher<geometry_msgs::msg::Point,SL_Bus_geometry_msgs_Point> Pub_nmpcRect_271;
// nmpcRect/NMPC Controller/State Publisher/Publish
SimulinkPublisher<catarob_interfaces::msg::StateEstimate,SL_Bus_catarob_interfaces_StateEstimate> Pub_nmpcRect_253;
// nmpcRect/Subscribe
SimulinkSubscriber<sensor_msgs::msg::NavSatFix,SL_Bus_sensor_msgs_NavSatFix> Sub_nmpcRect_272;
// nmpcRect/Subscribe1
SimulinkSubscriber<std_msgs::msg::Float64,SL_Bus_std_msgs_Float64> Sub_nmpcRect_273;
// nmpcRect/Subscribe2
SimulinkSubscriber<sensor_msgs::msg::Imu,SL_Bus_sensor_msgs_Imu> Sub_nmpcRect_274;
// For Block nmpcRect/Get Parameter
SimulinkParameterArrayGetter<real64_T,std::vector<double>> ParamGet_nmpcRect_289;

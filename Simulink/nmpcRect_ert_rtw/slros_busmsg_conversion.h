#ifndef _SLROS_BUSMSG_CONVERSION_H_
#define _SLROS_BUSMSG_CONVERSION_H_

#include "rclcpp/rclcpp.hpp"
#include <builtin_interfaces/msg/time.hpp>
#include <catarob_interfaces/msg/state_estimate.hpp>
#include <geometry_msgs/msg/point.hpp>
#include <geometry_msgs/msg/quaternion.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <geometry_msgs/msg/vector3.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/nav_sat_fix.hpp>
#include <sensor_msgs/msg/nav_sat_status.hpp>
#include <std_msgs/msg/float64.hpp>
#include <std_msgs/msg/header.hpp>
#include <std_msgs/msg/int32.hpp>
#include "nmpcRect_types.h"
#include "slros_msgconvert_utils.h"


void convertFromBus(builtin_interfaces::msg::Time& msgPtr, SL_Bus_builtin_interfaces_Time const* busPtr);
void convertToBus(SL_Bus_builtin_interfaces_Time* busPtr, const builtin_interfaces::msg::Time& msgPtr);

void convertFromBus(catarob_interfaces::msg::StateEstimate& msgPtr, SL_Bus_catarob_interfaces_StateEstimate const* busPtr);
void convertToBus(SL_Bus_catarob_interfaces_StateEstimate* busPtr, const catarob_interfaces::msg::StateEstimate& msgPtr);

void convertFromBus(geometry_msgs::msg::Point& msgPtr, SL_Bus_geometry_msgs_Point const* busPtr);
void convertToBus(SL_Bus_geometry_msgs_Point* busPtr, const geometry_msgs::msg::Point& msgPtr);

void convertFromBus(geometry_msgs::msg::Quaternion& msgPtr, SL_Bus_geometry_msgs_Quaternion const* busPtr);
void convertToBus(SL_Bus_geometry_msgs_Quaternion* busPtr, const geometry_msgs::msg::Quaternion& msgPtr);

void convertFromBus(geometry_msgs::msg::Twist& msgPtr, SL_Bus_geometry_msgs_Twist const* busPtr);
void convertToBus(SL_Bus_geometry_msgs_Twist* busPtr, const geometry_msgs::msg::Twist& msgPtr);

void convertFromBus(geometry_msgs::msg::Vector3& msgPtr, SL_Bus_geometry_msgs_Vector3 const* busPtr);
void convertToBus(SL_Bus_geometry_msgs_Vector3* busPtr, const geometry_msgs::msg::Vector3& msgPtr);

void convertFromBus(sensor_msgs::msg::Imu& msgPtr, SL_Bus_sensor_msgs_Imu const* busPtr);
void convertToBus(SL_Bus_sensor_msgs_Imu* busPtr, const sensor_msgs::msg::Imu& msgPtr);

void convertFromBus(sensor_msgs::msg::NavSatFix& msgPtr, SL_Bus_sensor_msgs_NavSatFix const* busPtr);
void convertToBus(SL_Bus_sensor_msgs_NavSatFix* busPtr, const sensor_msgs::msg::NavSatFix& msgPtr);

void convertFromBus(sensor_msgs::msg::NavSatStatus& msgPtr, SL_Bus_sensor_msgs_NavSatStatus const* busPtr);
void convertToBus(SL_Bus_sensor_msgs_NavSatStatus* busPtr, const sensor_msgs::msg::NavSatStatus& msgPtr);

void convertFromBus(std_msgs::msg::Float64& msgPtr, SL_Bus_std_msgs_Float64 const* busPtr);
void convertToBus(SL_Bus_std_msgs_Float64* busPtr, const std_msgs::msg::Float64& msgPtr);

void convertFromBus(std_msgs::msg::Header& msgPtr, SL_Bus_std_msgs_Header const* busPtr);
void convertToBus(SL_Bus_std_msgs_Header* busPtr, const std_msgs::msg::Header& msgPtr);

void convertFromBus(std_msgs::msg::Int32& msgPtr, SL_Bus_std_msgs_Int32 const* busPtr);
void convertToBus(SL_Bus_std_msgs_Int32* busPtr, const std_msgs::msg::Int32& msgPtr);


#endif

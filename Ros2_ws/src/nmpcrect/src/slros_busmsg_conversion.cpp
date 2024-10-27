#include "slros_busmsg_conversion.h"


// Conversions between SL_Bus_builtin_interfaces_Time and builtin_interfaces::msg::Time

void convertFromBus(builtin_interfaces::msg::Time& msgPtr, SL_Bus_builtin_interfaces_Time const* busPtr)
{
  const std::string rosMessageType("builtin_interfaces/Time");

  msgPtr.nanosec =  busPtr->nanosec;
  msgPtr.sec =  busPtr->sec;
}

void convertToBus(SL_Bus_builtin_interfaces_Time* busPtr, const builtin_interfaces::msg::Time& msgPtr)
{
  const std::string rosMessageType("builtin_interfaces/Time");

  busPtr->nanosec =  msgPtr.nanosec;
  busPtr->sec =  msgPtr.sec;
}


// Conversions between SL_Bus_catarob_interfaces_StateEstimate and catarob_interfaces::msg::StateEstimate

void convertFromBus(catarob_interfaces::msg::StateEstimate& msgPtr, SL_Bus_catarob_interfaces_StateEstimate const* busPtr)
{
  const std::string rosMessageType("catarob_interfaces/StateEstimate");

  convertFromBusFixedPrimitiveArray(msgPtr.covariance, busPtr->covariance);
  msgPtr.psi =  busPtr->psi;
  msgPtr.r =  busPtr->r;
  msgPtr.u =  busPtr->u;
  msgPtr.v =  busPtr->v;
  msgPtr.x =  busPtr->x;
  msgPtr.y =  busPtr->y;
}

void convertToBus(SL_Bus_catarob_interfaces_StateEstimate* busPtr, const catarob_interfaces::msg::StateEstimate& msgPtr)
{
  const std::string rosMessageType("catarob_interfaces/StateEstimate");

  convertToBusFixedPrimitiveArray(busPtr->covariance, msgPtr.covariance, slros::NoopWarning());
  busPtr->psi =  msgPtr.psi;
  busPtr->r =  msgPtr.r;
  busPtr->u =  msgPtr.u;
  busPtr->v =  msgPtr.v;
  busPtr->x =  msgPtr.x;
  busPtr->y =  msgPtr.y;
}


// Conversions between SL_Bus_geometry_msgs_Point and geometry_msgs::msg::Point

void convertFromBus(geometry_msgs::msg::Point& msgPtr, SL_Bus_geometry_msgs_Point const* busPtr)
{
  const std::string rosMessageType("geometry_msgs/Point");

  msgPtr.x =  busPtr->x;
  msgPtr.y =  busPtr->y;
  msgPtr.z =  busPtr->z;
}

void convertToBus(SL_Bus_geometry_msgs_Point* busPtr, const geometry_msgs::msg::Point& msgPtr)
{
  const std::string rosMessageType("geometry_msgs/Point");

  busPtr->x =  msgPtr.x;
  busPtr->y =  msgPtr.y;
  busPtr->z =  msgPtr.z;
}


// Conversions between SL_Bus_geometry_msgs_Quaternion and geometry_msgs::msg::Quaternion

void convertFromBus(geometry_msgs::msg::Quaternion& msgPtr, SL_Bus_geometry_msgs_Quaternion const* busPtr)
{
  const std::string rosMessageType("geometry_msgs/Quaternion");

  msgPtr.w =  busPtr->w;
  msgPtr.x =  busPtr->x;
  msgPtr.y =  busPtr->y;
  msgPtr.z =  busPtr->z;
}

void convertToBus(SL_Bus_geometry_msgs_Quaternion* busPtr, const geometry_msgs::msg::Quaternion& msgPtr)
{
  const std::string rosMessageType("geometry_msgs/Quaternion");

  busPtr->w =  msgPtr.w;
  busPtr->x =  msgPtr.x;
  busPtr->y =  msgPtr.y;
  busPtr->z =  msgPtr.z;
}


// Conversions between SL_Bus_geometry_msgs_Twist and geometry_msgs::msg::Twist

void convertFromBus(geometry_msgs::msg::Twist& msgPtr, SL_Bus_geometry_msgs_Twist const* busPtr)
{
  const std::string rosMessageType("geometry_msgs/Twist");

  convertFromBus(msgPtr.angular, &busPtr->angular);
  convertFromBus(msgPtr.linear, &busPtr->linear);
}

void convertToBus(SL_Bus_geometry_msgs_Twist* busPtr, const geometry_msgs::msg::Twist& msgPtr)
{
  const std::string rosMessageType("geometry_msgs/Twist");

  convertToBus(&busPtr->angular, msgPtr.angular);
  convertToBus(&busPtr->linear, msgPtr.linear);
}


// Conversions between SL_Bus_geometry_msgs_Vector3 and geometry_msgs::msg::Vector3

void convertFromBus(geometry_msgs::msg::Vector3& msgPtr, SL_Bus_geometry_msgs_Vector3 const* busPtr)
{
  const std::string rosMessageType("geometry_msgs/Vector3");

  msgPtr.x =  busPtr->x;
  msgPtr.y =  busPtr->y;
  msgPtr.z =  busPtr->z;
}

void convertToBus(SL_Bus_geometry_msgs_Vector3* busPtr, const geometry_msgs::msg::Vector3& msgPtr)
{
  const std::string rosMessageType("geometry_msgs/Vector3");

  busPtr->x =  msgPtr.x;
  busPtr->y =  msgPtr.y;
  busPtr->z =  msgPtr.z;
}


// Conversions between SL_Bus_sensor_msgs_Imu and sensor_msgs::msg::Imu

void convertFromBus(sensor_msgs::msg::Imu& msgPtr, SL_Bus_sensor_msgs_Imu const* busPtr)
{
  const std::string rosMessageType("sensor_msgs/Imu");

  convertFromBus(msgPtr.angular_velocity, &busPtr->angular_velocity);
  convertFromBusFixedPrimitiveArray(msgPtr.angular_velocity_covariance, busPtr->angular_velocity_covariance);
  convertFromBus(msgPtr.header, &busPtr->header);
  convertFromBus(msgPtr.linear_acceleration, &busPtr->linear_acceleration);
  convertFromBusFixedPrimitiveArray(msgPtr.linear_acceleration_covariance, busPtr->linear_acceleration_covariance);
  convertFromBus(msgPtr.orientation, &busPtr->orientation);
  convertFromBusFixedPrimitiveArray(msgPtr.orientation_covariance, busPtr->orientation_covariance);
}

void convertToBus(SL_Bus_sensor_msgs_Imu* busPtr, const sensor_msgs::msg::Imu& msgPtr)
{
  const std::string rosMessageType("sensor_msgs/Imu");

  convertToBus(&busPtr->angular_velocity, msgPtr.angular_velocity);
  convertToBusFixedPrimitiveArray(busPtr->angular_velocity_covariance, msgPtr.angular_velocity_covariance, slros::NoopWarning());
  convertToBus(&busPtr->header, msgPtr.header);
  convertToBus(&busPtr->linear_acceleration, msgPtr.linear_acceleration);
  convertToBusFixedPrimitiveArray(busPtr->linear_acceleration_covariance, msgPtr.linear_acceleration_covariance, slros::NoopWarning());
  convertToBus(&busPtr->orientation, msgPtr.orientation);
  convertToBusFixedPrimitiveArray(busPtr->orientation_covariance, msgPtr.orientation_covariance, slros::NoopWarning());
}


// Conversions between SL_Bus_sensor_msgs_NavSatFix and sensor_msgs::msg::NavSatFix

void convertFromBus(sensor_msgs::msg::NavSatFix& msgPtr, SL_Bus_sensor_msgs_NavSatFix const* busPtr)
{
  const std::string rosMessageType("sensor_msgs/NavSatFix");

  msgPtr.altitude =  busPtr->altitude;
  convertFromBus(msgPtr.header, &busPtr->header);
  msgPtr.latitude =  busPtr->latitude;
  msgPtr.longitude =  busPtr->longitude;
  convertFromBusFixedPrimitiveArray(msgPtr.position_covariance, busPtr->position_covariance);
  msgPtr.position_covariance_type =  busPtr->position_covariance_type;
  convertFromBus(msgPtr.status, &busPtr->status);
}

void convertToBus(SL_Bus_sensor_msgs_NavSatFix* busPtr, const sensor_msgs::msg::NavSatFix& msgPtr)
{
  const std::string rosMessageType("sensor_msgs/NavSatFix");

  busPtr->altitude =  msgPtr.altitude;
  convertToBus(&busPtr->header, msgPtr.header);
  busPtr->latitude =  msgPtr.latitude;
  busPtr->longitude =  msgPtr.longitude;
  convertToBusFixedPrimitiveArray(busPtr->position_covariance, msgPtr.position_covariance, slros::NoopWarning());
  busPtr->position_covariance_type =  msgPtr.position_covariance_type;
  convertToBus(&busPtr->status, msgPtr.status);
}


// Conversions between SL_Bus_sensor_msgs_NavSatStatus and sensor_msgs::msg::NavSatStatus

void convertFromBus(sensor_msgs::msg::NavSatStatus& msgPtr, SL_Bus_sensor_msgs_NavSatStatus const* busPtr)
{
  const std::string rosMessageType("sensor_msgs/NavSatStatus");

  msgPtr.service =  busPtr->service;
  msgPtr.status =  busPtr->status;
}

void convertToBus(SL_Bus_sensor_msgs_NavSatStatus* busPtr, const sensor_msgs::msg::NavSatStatus& msgPtr)
{
  const std::string rosMessageType("sensor_msgs/NavSatStatus");

  busPtr->service =  msgPtr.service;
  busPtr->status =  msgPtr.status;
}


// Conversions between SL_Bus_std_msgs_Float64 and std_msgs::msg::Float64

void convertFromBus(std_msgs::msg::Float64& msgPtr, SL_Bus_std_msgs_Float64 const* busPtr)
{
  const std::string rosMessageType("std_msgs/Float64");

  msgPtr.data =  busPtr->data;
}

void convertToBus(SL_Bus_std_msgs_Float64* busPtr, const std_msgs::msg::Float64& msgPtr)
{
  const std::string rosMessageType("std_msgs/Float64");

  busPtr->data =  msgPtr.data;
}


// Conversions between SL_Bus_std_msgs_Header and std_msgs::msg::Header

void convertFromBus(std_msgs::msg::Header& msgPtr, SL_Bus_std_msgs_Header const* busPtr)
{
  const std::string rosMessageType("std_msgs/Header");

  convertFromBusVariablePrimitiveArray(msgPtr.frame_id, busPtr->frame_id, busPtr->frame_id_SL_Info);
  convertFromBus(msgPtr.stamp, &busPtr->stamp);
}

void convertToBus(SL_Bus_std_msgs_Header* busPtr, const std_msgs::msg::Header& msgPtr)
{
  const std::string rosMessageType("std_msgs/Header");

  convertToBusVariablePrimitiveArray(busPtr->frame_id, busPtr->frame_id_SL_Info, msgPtr.frame_id, slros::EnabledWarning(rosMessageType, "frame_id"));
  convertToBus(&busPtr->stamp, msgPtr.stamp);
}


// Conversions between SL_Bus_std_msgs_Int32 and std_msgs::msg::Int32

void convertFromBus(std_msgs::msg::Int32& msgPtr, SL_Bus_std_msgs_Int32 const* busPtr)
{
  const std::string rosMessageType("std_msgs/Int32");

  msgPtr.data =  busPtr->data;
}

void convertToBus(SL_Bus_std_msgs_Int32* busPtr, const std_msgs::msg::Int32& msgPtr)
{
  const std::string rosMessageType("std_msgs/Int32");

  busPtr->data =  msgPtr.data;
}


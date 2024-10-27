// generated from rosidl_generator_cpp/resource/idl__builder.hpp.em
// with input from catarob_interfaces:msg\HullStatus.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__BUILDER_HPP_
#define CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__BUILDER_HPP_

#include <algorithm>
#include <utility>

#include "catarob_interfaces/msg/detail/hull_status__struct.hpp"
#include "rosidl_runtime_cpp/message_initialization.hpp"


namespace catarob_interfaces
{

namespace msg
{

namespace builder
{

class Init_HullStatus_pwr_relay
{
public:
  explicit Init_HullStatus_pwr_relay(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  ::catarob_interfaces::msg::HullStatus pwr_relay(::catarob_interfaces::msg::HullStatus::_pwr_relay_type arg)
  {
    msg_.pwr_relay = std::move(arg);
    return std::move(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_water_ingress
{
public:
  explicit Init_HullStatus_water_ingress(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_pwr_relay water_ingress(::catarob_interfaces::msg::HullStatus::_water_ingress_type arg)
  {
    msg_.water_ingress = std::move(arg);
    return Init_HullStatus_pwr_relay(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_tor_rc
{
public:
  explicit Init_HullStatus_tor_rc(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_water_ingress tor_rc(::catarob_interfaces::msg::HullStatus::_tor_rc_type arg)
  {
    msg_.tor_rc = std::move(arg);
    return Init_HullStatus_water_ingress(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_actual_imax
{
public:
  explicit Init_HullStatus_actual_imax(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_tor_rc actual_imax(::catarob_interfaces::msg::HullStatus::_actual_imax_type arg)
  {
    msg_.actual_imax = std::move(arg);
    return Init_HullStatus_tor_rc(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_vbat
{
public:
  explicit Init_HullStatus_vbat(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_actual_imax vbat(::catarob_interfaces::msg::HullStatus::_vbat_type arg)
  {
    msg_.vbat = std::move(arg);
    return Init_HullStatus_actual_imax(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_imot
{
public:
  explicit Init_HullStatus_imot(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_vbat imot(::catarob_interfaces::msg::HullStatus::_imot_type arg)
  {
    msg_.imot = std::move(arg);
    return Init_HullStatus_vbat(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_ibat
{
public:
  explicit Init_HullStatus_ibat(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_imot ibat(::catarob_interfaces::msg::HullStatus::_ibat_type arg)
  {
    msg_.ibat = std::move(arg);
    return Init_HullStatus_imot(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_pwm_source
{
public:
  explicit Init_HullStatus_pwm_source(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_ibat pwm_source(::catarob_interfaces::msg::HullStatus::_pwm_source_type arg)
  {
    msg_.pwm_source = std::move(arg);
    return Init_HullStatus_ibat(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_actual_pwm_rc
{
public:
  explicit Init_HullStatus_actual_pwm_rc(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_pwm_source actual_pwm_rc(::catarob_interfaces::msg::HullStatus::_actual_pwm_rc_type arg)
  {
    msg_.actual_pwm_rc = std::move(arg);
    return Init_HullStatus_pwm_source(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_actual_pwm
{
public:
  explicit Init_HullStatus_actual_pwm(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_actual_pwm_rc actual_pwm(::catarob_interfaces::msg::HullStatus::_actual_pwm_type arg)
  {
    msg_.actual_pwm = std::move(arg);
    return Init_HullStatus_actual_pwm_rc(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_temp
{
public:
  explicit Init_HullStatus_temp(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_actual_pwm temp(::catarob_interfaces::msg::HullStatus::_temp_type arg)
  {
    msg_.temp = std::move(arg);
    return Init_HullStatus_actual_pwm(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_status
{
public:
  explicit Init_HullStatus_status(::catarob_interfaces::msg::HullStatus & msg)
  : msg_(msg)
  {}
  Init_HullStatus_temp status(::catarob_interfaces::msg::HullStatus::_status_type arg)
  {
    msg_.status = std::move(arg);
    return Init_HullStatus_temp(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

class Init_HullStatus_header
{
public:
  Init_HullStatus_header()
  : msg_(::rosidl_runtime_cpp::MessageInitialization::SKIP)
  {}
  Init_HullStatus_status header(::catarob_interfaces::msg::HullStatus::_header_type arg)
  {
    msg_.header = std::move(arg);
    return Init_HullStatus_status(msg_);
  }

private:
  ::catarob_interfaces::msg::HullStatus msg_;
};

}  // namespace builder

}  // namespace msg

template<typename MessageType>
auto build();

template<>
inline
auto build<::catarob_interfaces::msg::HullStatus>()
{
  return catarob_interfaces::msg::builder::Init_HullStatus_header();
}

}  // namespace catarob_interfaces

#endif  // CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__BUILDER_HPP_

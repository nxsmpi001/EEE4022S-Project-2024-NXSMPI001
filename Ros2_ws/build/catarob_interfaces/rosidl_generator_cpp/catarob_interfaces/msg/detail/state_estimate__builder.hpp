// generated from rosidl_generator_cpp/resource/idl__builder.hpp.em
// with input from catarob_interfaces:msg/StateEstimate.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__BUILDER_HPP_
#define CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__BUILDER_HPP_

#include <algorithm>
#include <utility>

#include "catarob_interfaces/msg/detail/state_estimate__struct.hpp"
#include "rosidl_runtime_cpp/message_initialization.hpp"


namespace catarob_interfaces
{

namespace msg
{

namespace builder
{

class Init_StateEstimate_covariance
{
public:
  explicit Init_StateEstimate_covariance(::catarob_interfaces::msg::StateEstimate & msg)
  : msg_(msg)
  {}
  ::catarob_interfaces::msg::StateEstimate covariance(::catarob_interfaces::msg::StateEstimate::_covariance_type arg)
  {
    msg_.covariance = std::move(arg);
    return std::move(msg_);
  }

private:
  ::catarob_interfaces::msg::StateEstimate msg_;
};

class Init_StateEstimate_r
{
public:
  explicit Init_StateEstimate_r(::catarob_interfaces::msg::StateEstimate & msg)
  : msg_(msg)
  {}
  Init_StateEstimate_covariance r(::catarob_interfaces::msg::StateEstimate::_r_type arg)
  {
    msg_.r = std::move(arg);
    return Init_StateEstimate_covariance(msg_);
  }

private:
  ::catarob_interfaces::msg::StateEstimate msg_;
};

class Init_StateEstimate_v
{
public:
  explicit Init_StateEstimate_v(::catarob_interfaces::msg::StateEstimate & msg)
  : msg_(msg)
  {}
  Init_StateEstimate_r v(::catarob_interfaces::msg::StateEstimate::_v_type arg)
  {
    msg_.v = std::move(arg);
    return Init_StateEstimate_r(msg_);
  }

private:
  ::catarob_interfaces::msg::StateEstimate msg_;
};

class Init_StateEstimate_u
{
public:
  explicit Init_StateEstimate_u(::catarob_interfaces::msg::StateEstimate & msg)
  : msg_(msg)
  {}
  Init_StateEstimate_v u(::catarob_interfaces::msg::StateEstimate::_u_type arg)
  {
    msg_.u = std::move(arg);
    return Init_StateEstimate_v(msg_);
  }

private:
  ::catarob_interfaces::msg::StateEstimate msg_;
};

class Init_StateEstimate_psi
{
public:
  explicit Init_StateEstimate_psi(::catarob_interfaces::msg::StateEstimate & msg)
  : msg_(msg)
  {}
  Init_StateEstimate_u psi(::catarob_interfaces::msg::StateEstimate::_psi_type arg)
  {
    msg_.psi = std::move(arg);
    return Init_StateEstimate_u(msg_);
  }

private:
  ::catarob_interfaces::msg::StateEstimate msg_;
};

class Init_StateEstimate_y
{
public:
  explicit Init_StateEstimate_y(::catarob_interfaces::msg::StateEstimate & msg)
  : msg_(msg)
  {}
  Init_StateEstimate_psi y(::catarob_interfaces::msg::StateEstimate::_y_type arg)
  {
    msg_.y = std::move(arg);
    return Init_StateEstimate_psi(msg_);
  }

private:
  ::catarob_interfaces::msg::StateEstimate msg_;
};

class Init_StateEstimate_x
{
public:
  Init_StateEstimate_x()
  : msg_(::rosidl_runtime_cpp::MessageInitialization::SKIP)
  {}
  Init_StateEstimate_y x(::catarob_interfaces::msg::StateEstimate::_x_type arg)
  {
    msg_.x = std::move(arg);
    return Init_StateEstimate_y(msg_);
  }

private:
  ::catarob_interfaces::msg::StateEstimate msg_;
};

}  // namespace builder

}  // namespace msg

template<typename MessageType>
auto build();

template<>
inline
auto build<::catarob_interfaces::msg::StateEstimate>()
{
  return catarob_interfaces::msg::builder::Init_StateEstimate_x();
}

}  // namespace catarob_interfaces

#endif  // CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__BUILDER_HPP_

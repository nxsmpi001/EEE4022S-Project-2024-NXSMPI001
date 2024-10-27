// generated from rosidl_generator_cpp/resource/idl__builder.hpp.em
// with input from catarob_interfaces:srv\ToggleProps.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__BUILDER_HPP_
#define CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__BUILDER_HPP_

#include <algorithm>
#include <utility>

#include "catarob_interfaces/srv/detail/toggle_props__struct.hpp"
#include "rosidl_runtime_cpp/message_initialization.hpp"


namespace catarob_interfaces
{

namespace srv
{

namespace builder
{

class Init_ToggleProps_Request_flag
{
public:
  Init_ToggleProps_Request_flag()
  : msg_(::rosidl_runtime_cpp::MessageInitialization::SKIP)
  {}
  ::catarob_interfaces::srv::ToggleProps_Request flag(::catarob_interfaces::srv::ToggleProps_Request::_flag_type arg)
  {
    msg_.flag = std::move(arg);
    return std::move(msg_);
  }

private:
  ::catarob_interfaces::srv::ToggleProps_Request msg_;
};

}  // namespace builder

}  // namespace srv

template<typename MessageType>
auto build();

template<>
inline
auto build<::catarob_interfaces::srv::ToggleProps_Request>()
{
  return catarob_interfaces::srv::builder::Init_ToggleProps_Request_flag();
}

}  // namespace catarob_interfaces


namespace catarob_interfaces
{

namespace srv
{

namespace builder
{

class Init_ToggleProps_Response_state
{
public:
  Init_ToggleProps_Response_state()
  : msg_(::rosidl_runtime_cpp::MessageInitialization::SKIP)
  {}
  ::catarob_interfaces::srv::ToggleProps_Response state(::catarob_interfaces::srv::ToggleProps_Response::_state_type arg)
  {
    msg_.state = std::move(arg);
    return std::move(msg_);
  }

private:
  ::catarob_interfaces::srv::ToggleProps_Response msg_;
};

}  // namespace builder

}  // namespace srv

template<typename MessageType>
auto build();

template<>
inline
auto build<::catarob_interfaces::srv::ToggleProps_Response>()
{
  return catarob_interfaces::srv::builder::Init_ToggleProps_Response_state();
}

}  // namespace catarob_interfaces

#endif  // CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__BUILDER_HPP_

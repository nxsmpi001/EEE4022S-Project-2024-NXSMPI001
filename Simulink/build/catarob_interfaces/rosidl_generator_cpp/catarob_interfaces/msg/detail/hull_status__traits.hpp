// generated from rosidl_generator_cpp/resource/idl__traits.hpp.em
// with input from catarob_interfaces:msg\HullStatus.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__TRAITS_HPP_
#define CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__TRAITS_HPP_

#include <stdint.h>

#include <sstream>
#include <string>
#include <type_traits>

#include "catarob_interfaces/msg/detail/hull_status__struct.hpp"
#include "rosidl_runtime_cpp/traits.hpp"

// Include directives for member types
// Member 'header'
#include "std_msgs/msg/detail/header__traits.hpp"

namespace catarob_interfaces
{

namespace msg
{

inline void to_flow_style_yaml(
  const HullStatus & msg,
  std::ostream & out)
{
  out << "{";
  // member: header
  {
    out << "header: ";
    to_flow_style_yaml(msg.header, out);
    out << ", ";
  }

  // member: status
  {
    out << "status: ";
    rosidl_generator_traits::value_to_yaml(msg.status, out);
    out << ", ";
  }

  // member: temp
  {
    out << "temp: ";
    rosidl_generator_traits::value_to_yaml(msg.temp, out);
    out << ", ";
  }

  // member: actual_pwm
  {
    out << "actual_pwm: ";
    rosidl_generator_traits::value_to_yaml(msg.actual_pwm, out);
    out << ", ";
  }

  // member: actual_pwm_rc
  {
    out << "actual_pwm_rc: ";
    rosidl_generator_traits::value_to_yaml(msg.actual_pwm_rc, out);
    out << ", ";
  }

  // member: pwm_source
  {
    out << "pwm_source: ";
    rosidl_generator_traits::value_to_yaml(msg.pwm_source, out);
    out << ", ";
  }

  // member: ibat
  {
    out << "ibat: ";
    rosidl_generator_traits::value_to_yaml(msg.ibat, out);
    out << ", ";
  }

  // member: imot
  {
    out << "imot: ";
    rosidl_generator_traits::value_to_yaml(msg.imot, out);
    out << ", ";
  }

  // member: vbat
  {
    out << "vbat: ";
    rosidl_generator_traits::value_to_yaml(msg.vbat, out);
    out << ", ";
  }

  // member: actual_imax
  {
    out << "actual_imax: ";
    rosidl_generator_traits::value_to_yaml(msg.actual_imax, out);
    out << ", ";
  }

  // member: tor_rc
  {
    out << "tor_rc: ";
    rosidl_generator_traits::value_to_yaml(msg.tor_rc, out);
    out << ", ";
  }

  // member: water_ingress
  {
    out << "water_ingress: ";
    rosidl_generator_traits::value_to_yaml(msg.water_ingress, out);
    out << ", ";
  }

  // member: pwr_relay
  {
    out << "pwr_relay: ";
    rosidl_generator_traits::value_to_yaml(msg.pwr_relay, out);
  }
  out << "}";
}  // NOLINT(readability/fn_size)

inline void to_block_style_yaml(
  const HullStatus & msg,
  std::ostream & out, size_t indentation = 0)
{
  // member: header
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "header:\n";
    to_block_style_yaml(msg.header, out, indentation + 2);
  }

  // member: status
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "status: ";
    rosidl_generator_traits::value_to_yaml(msg.status, out);
    out << "\n";
  }

  // member: temp
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "temp: ";
    rosidl_generator_traits::value_to_yaml(msg.temp, out);
    out << "\n";
  }

  // member: actual_pwm
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "actual_pwm: ";
    rosidl_generator_traits::value_to_yaml(msg.actual_pwm, out);
    out << "\n";
  }

  // member: actual_pwm_rc
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "actual_pwm_rc: ";
    rosidl_generator_traits::value_to_yaml(msg.actual_pwm_rc, out);
    out << "\n";
  }

  // member: pwm_source
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "pwm_source: ";
    rosidl_generator_traits::value_to_yaml(msg.pwm_source, out);
    out << "\n";
  }

  // member: ibat
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "ibat: ";
    rosidl_generator_traits::value_to_yaml(msg.ibat, out);
    out << "\n";
  }

  // member: imot
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "imot: ";
    rosidl_generator_traits::value_to_yaml(msg.imot, out);
    out << "\n";
  }

  // member: vbat
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "vbat: ";
    rosidl_generator_traits::value_to_yaml(msg.vbat, out);
    out << "\n";
  }

  // member: actual_imax
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "actual_imax: ";
    rosidl_generator_traits::value_to_yaml(msg.actual_imax, out);
    out << "\n";
  }

  // member: tor_rc
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "tor_rc: ";
    rosidl_generator_traits::value_to_yaml(msg.tor_rc, out);
    out << "\n";
  }

  // member: water_ingress
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "water_ingress: ";
    rosidl_generator_traits::value_to_yaml(msg.water_ingress, out);
    out << "\n";
  }

  // member: pwr_relay
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "pwr_relay: ";
    rosidl_generator_traits::value_to_yaml(msg.pwr_relay, out);
    out << "\n";
  }
}  // NOLINT(readability/fn_size)

inline std::string to_yaml(const HullStatus & msg, bool use_flow_style = false)
{
  std::ostringstream out;
  if (use_flow_style) {
    to_flow_style_yaml(msg, out);
  } else {
    to_block_style_yaml(msg, out);
  }
  return out.str();
}

}  // namespace msg

}  // namespace catarob_interfaces

namespace rosidl_generator_traits
{

[[deprecated("use catarob_interfaces::msg::to_block_style_yaml() instead")]]
inline void to_yaml(
  const catarob_interfaces::msg::HullStatus & msg,
  std::ostream & out, size_t indentation = 0)
{
  catarob_interfaces::msg::to_block_style_yaml(msg, out, indentation);
}

[[deprecated("use catarob_interfaces::msg::to_yaml() instead")]]
inline std::string to_yaml(const catarob_interfaces::msg::HullStatus & msg)
{
  return catarob_interfaces::msg::to_yaml(msg);
}

template<>
inline const char * data_type<catarob_interfaces::msg::HullStatus>()
{
  return "catarob_interfaces::msg::HullStatus";
}

template<>
inline const char * name<catarob_interfaces::msg::HullStatus>()
{
  return "catarob_interfaces/msg/HullStatus";
}

template<>
struct has_fixed_size<catarob_interfaces::msg::HullStatus>
  : std::integral_constant<bool, has_fixed_size<std_msgs::msg::Header>::value> {};

template<>
struct has_bounded_size<catarob_interfaces::msg::HullStatus>
  : std::integral_constant<bool, has_bounded_size<std_msgs::msg::Header>::value> {};

template<>
struct is_message<catarob_interfaces::msg::HullStatus>
  : std::true_type {};

}  // namespace rosidl_generator_traits

#endif  // CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__TRAITS_HPP_

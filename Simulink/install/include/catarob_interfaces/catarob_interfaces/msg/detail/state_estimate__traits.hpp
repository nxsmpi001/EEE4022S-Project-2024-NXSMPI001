// generated from rosidl_generator_cpp/resource/idl__traits.hpp.em
// with input from catarob_interfaces:msg\StateEstimate.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__TRAITS_HPP_
#define CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__TRAITS_HPP_

#include <stdint.h>

#include <sstream>
#include <string>
#include <type_traits>

#include "catarob_interfaces/msg/detail/state_estimate__struct.hpp"
#include "rosidl_runtime_cpp/traits.hpp"

namespace catarob_interfaces
{

namespace msg
{

inline void to_flow_style_yaml(
  const StateEstimate & msg,
  std::ostream & out)
{
  out << "{";
  // member: x
  {
    out << "x: ";
    rosidl_generator_traits::value_to_yaml(msg.x, out);
    out << ", ";
  }

  // member: y
  {
    out << "y: ";
    rosidl_generator_traits::value_to_yaml(msg.y, out);
    out << ", ";
  }

  // member: psi
  {
    out << "psi: ";
    rosidl_generator_traits::value_to_yaml(msg.psi, out);
    out << ", ";
  }

  // member: u
  {
    out << "u: ";
    rosidl_generator_traits::value_to_yaml(msg.u, out);
    out << ", ";
  }

  // member: v
  {
    out << "v: ";
    rosidl_generator_traits::value_to_yaml(msg.v, out);
    out << ", ";
  }

  // member: r
  {
    out << "r: ";
    rosidl_generator_traits::value_to_yaml(msg.r, out);
    out << ", ";
  }

  // member: covariance
  {
    if (msg.covariance.size() == 0) {
      out << "covariance: []";
    } else {
      out << "covariance: [";
      size_t pending_items = msg.covariance.size();
      for (auto item : msg.covariance) {
        rosidl_generator_traits::value_to_yaml(item, out);
        if (--pending_items > 0) {
          out << ", ";
        }
      }
      out << "]";
    }
  }
  out << "}";
}  // NOLINT(readability/fn_size)

inline void to_block_style_yaml(
  const StateEstimate & msg,
  std::ostream & out, size_t indentation = 0)
{
  // member: x
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "x: ";
    rosidl_generator_traits::value_to_yaml(msg.x, out);
    out << "\n";
  }

  // member: y
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "y: ";
    rosidl_generator_traits::value_to_yaml(msg.y, out);
    out << "\n";
  }

  // member: psi
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "psi: ";
    rosidl_generator_traits::value_to_yaml(msg.psi, out);
    out << "\n";
  }

  // member: u
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "u: ";
    rosidl_generator_traits::value_to_yaml(msg.u, out);
    out << "\n";
  }

  // member: v
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "v: ";
    rosidl_generator_traits::value_to_yaml(msg.v, out);
    out << "\n";
  }

  // member: r
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "r: ";
    rosidl_generator_traits::value_to_yaml(msg.r, out);
    out << "\n";
  }

  // member: covariance
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    if (msg.covariance.size() == 0) {
      out << "covariance: []\n";
    } else {
      out << "covariance:\n";
      for (auto item : msg.covariance) {
        if (indentation > 0) {
          out << std::string(indentation, ' ');
        }
        out << "- ";
        rosidl_generator_traits::value_to_yaml(item, out);
        out << "\n";
      }
    }
  }
}  // NOLINT(readability/fn_size)

inline std::string to_yaml(const StateEstimate & msg, bool use_flow_style = false)
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
  const catarob_interfaces::msg::StateEstimate & msg,
  std::ostream & out, size_t indentation = 0)
{
  catarob_interfaces::msg::to_block_style_yaml(msg, out, indentation);
}

[[deprecated("use catarob_interfaces::msg::to_yaml() instead")]]
inline std::string to_yaml(const catarob_interfaces::msg::StateEstimate & msg)
{
  return catarob_interfaces::msg::to_yaml(msg);
}

template<>
inline const char * data_type<catarob_interfaces::msg::StateEstimate>()
{
  return "catarob_interfaces::msg::StateEstimate";
}

template<>
inline const char * name<catarob_interfaces::msg::StateEstimate>()
{
  return "catarob_interfaces/msg/StateEstimate";
}

template<>
struct has_fixed_size<catarob_interfaces::msg::StateEstimate>
  : std::integral_constant<bool, true> {};

template<>
struct has_bounded_size<catarob_interfaces::msg::StateEstimate>
  : std::integral_constant<bool, true> {};

template<>
struct is_message<catarob_interfaces::msg::StateEstimate>
  : std::true_type {};

}  // namespace rosidl_generator_traits

#endif  // CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__TRAITS_HPP_

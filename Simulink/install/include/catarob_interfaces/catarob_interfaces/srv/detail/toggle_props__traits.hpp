// generated from rosidl_generator_cpp/resource/idl__traits.hpp.em
// with input from catarob_interfaces:srv\ToggleProps.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__TRAITS_HPP_
#define CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__TRAITS_HPP_

#include <stdint.h>

#include <sstream>
#include <string>
#include <type_traits>

#include "catarob_interfaces/srv/detail/toggle_props__struct.hpp"
#include "rosidl_runtime_cpp/traits.hpp"

namespace catarob_interfaces
{

namespace srv
{

inline void to_flow_style_yaml(
  const ToggleProps_Request & msg,
  std::ostream & out)
{
  out << "{";
  // member: flag
  {
    out << "flag: ";
    rosidl_generator_traits::value_to_yaml(msg.flag, out);
  }
  out << "}";
}  // NOLINT(readability/fn_size)

inline void to_block_style_yaml(
  const ToggleProps_Request & msg,
  std::ostream & out, size_t indentation = 0)
{
  // member: flag
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "flag: ";
    rosidl_generator_traits::value_to_yaml(msg.flag, out);
    out << "\n";
  }
}  // NOLINT(readability/fn_size)

inline std::string to_yaml(const ToggleProps_Request & msg, bool use_flow_style = false)
{
  std::ostringstream out;
  if (use_flow_style) {
    to_flow_style_yaml(msg, out);
  } else {
    to_block_style_yaml(msg, out);
  }
  return out.str();
}

}  // namespace srv

}  // namespace catarob_interfaces

namespace rosidl_generator_traits
{

[[deprecated("use catarob_interfaces::srv::to_block_style_yaml() instead")]]
inline void to_yaml(
  const catarob_interfaces::srv::ToggleProps_Request & msg,
  std::ostream & out, size_t indentation = 0)
{
  catarob_interfaces::srv::to_block_style_yaml(msg, out, indentation);
}

[[deprecated("use catarob_interfaces::srv::to_yaml() instead")]]
inline std::string to_yaml(const catarob_interfaces::srv::ToggleProps_Request & msg)
{
  return catarob_interfaces::srv::to_yaml(msg);
}

template<>
inline const char * data_type<catarob_interfaces::srv::ToggleProps_Request>()
{
  return "catarob_interfaces::srv::ToggleProps_Request";
}

template<>
inline const char * name<catarob_interfaces::srv::ToggleProps_Request>()
{
  return "catarob_interfaces/srv/ToggleProps_Request";
}

template<>
struct has_fixed_size<catarob_interfaces::srv::ToggleProps_Request>
  : std::integral_constant<bool, true> {};

template<>
struct has_bounded_size<catarob_interfaces::srv::ToggleProps_Request>
  : std::integral_constant<bool, true> {};

template<>
struct is_message<catarob_interfaces::srv::ToggleProps_Request>
  : std::true_type {};

}  // namespace rosidl_generator_traits

namespace catarob_interfaces
{

namespace srv
{

inline void to_flow_style_yaml(
  const ToggleProps_Response & msg,
  std::ostream & out)
{
  out << "{";
  // member: state
  {
    out << "state: ";
    rosidl_generator_traits::value_to_yaml(msg.state, out);
  }
  out << "}";
}  // NOLINT(readability/fn_size)

inline void to_block_style_yaml(
  const ToggleProps_Response & msg,
  std::ostream & out, size_t indentation = 0)
{
  // member: state
  {
    if (indentation > 0) {
      out << std::string(indentation, ' ');
    }
    out << "state: ";
    rosidl_generator_traits::value_to_yaml(msg.state, out);
    out << "\n";
  }
}  // NOLINT(readability/fn_size)

inline std::string to_yaml(const ToggleProps_Response & msg, bool use_flow_style = false)
{
  std::ostringstream out;
  if (use_flow_style) {
    to_flow_style_yaml(msg, out);
  } else {
    to_block_style_yaml(msg, out);
  }
  return out.str();
}

}  // namespace srv

}  // namespace catarob_interfaces

namespace rosidl_generator_traits
{

[[deprecated("use catarob_interfaces::srv::to_block_style_yaml() instead")]]
inline void to_yaml(
  const catarob_interfaces::srv::ToggleProps_Response & msg,
  std::ostream & out, size_t indentation = 0)
{
  catarob_interfaces::srv::to_block_style_yaml(msg, out, indentation);
}

[[deprecated("use catarob_interfaces::srv::to_yaml() instead")]]
inline std::string to_yaml(const catarob_interfaces::srv::ToggleProps_Response & msg)
{
  return catarob_interfaces::srv::to_yaml(msg);
}

template<>
inline const char * data_type<catarob_interfaces::srv::ToggleProps_Response>()
{
  return "catarob_interfaces::srv::ToggleProps_Response";
}

template<>
inline const char * name<catarob_interfaces::srv::ToggleProps_Response>()
{
  return "catarob_interfaces/srv/ToggleProps_Response";
}

template<>
struct has_fixed_size<catarob_interfaces::srv::ToggleProps_Response>
  : std::integral_constant<bool, true> {};

template<>
struct has_bounded_size<catarob_interfaces::srv::ToggleProps_Response>
  : std::integral_constant<bool, true> {};

template<>
struct is_message<catarob_interfaces::srv::ToggleProps_Response>
  : std::true_type {};

}  // namespace rosidl_generator_traits

namespace rosidl_generator_traits
{

template<>
inline const char * data_type<catarob_interfaces::srv::ToggleProps>()
{
  return "catarob_interfaces::srv::ToggleProps";
}

template<>
inline const char * name<catarob_interfaces::srv::ToggleProps>()
{
  return "catarob_interfaces/srv/ToggleProps";
}

template<>
struct has_fixed_size<catarob_interfaces::srv::ToggleProps>
  : std::integral_constant<
    bool,
    has_fixed_size<catarob_interfaces::srv::ToggleProps_Request>::value &&
    has_fixed_size<catarob_interfaces::srv::ToggleProps_Response>::value
  >
{
};

template<>
struct has_bounded_size<catarob_interfaces::srv::ToggleProps>
  : std::integral_constant<
    bool,
    has_bounded_size<catarob_interfaces::srv::ToggleProps_Request>::value &&
    has_bounded_size<catarob_interfaces::srv::ToggleProps_Response>::value
  >
{
};

template<>
struct is_service<catarob_interfaces::srv::ToggleProps>
  : std::true_type
{
};

template<>
struct is_service_request<catarob_interfaces::srv::ToggleProps_Request>
  : std::true_type
{
};

template<>
struct is_service_response<catarob_interfaces::srv::ToggleProps_Response>
  : std::true_type
{
};

}  // namespace rosidl_generator_traits

#endif  // CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__TRAITS_HPP_

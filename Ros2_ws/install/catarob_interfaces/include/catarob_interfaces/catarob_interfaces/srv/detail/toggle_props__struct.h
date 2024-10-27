// generated from rosidl_generator_c/resource/idl__struct.h.em
// with input from catarob_interfaces:srv/ToggleProps.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__STRUCT_H_
#define CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__STRUCT_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>


// Constants defined in the message

/// Struct defined in srv/ToggleProps in the package catarob_interfaces.
typedef struct catarob_interfaces__srv__ToggleProps_Request
{
  bool flag;
} catarob_interfaces__srv__ToggleProps_Request;

// Struct for a sequence of catarob_interfaces__srv__ToggleProps_Request.
typedef struct catarob_interfaces__srv__ToggleProps_Request__Sequence
{
  catarob_interfaces__srv__ToggleProps_Request * data;
  /// The number of valid items in data
  size_t size;
  /// The number of allocated items in data
  size_t capacity;
} catarob_interfaces__srv__ToggleProps_Request__Sequence;


// Constants defined in the message

/// Struct defined in srv/ToggleProps in the package catarob_interfaces.
typedef struct catarob_interfaces__srv__ToggleProps_Response
{
  bool state;
} catarob_interfaces__srv__ToggleProps_Response;

// Struct for a sequence of catarob_interfaces__srv__ToggleProps_Response.
typedef struct catarob_interfaces__srv__ToggleProps_Response__Sequence
{
  catarob_interfaces__srv__ToggleProps_Response * data;
  /// The number of valid items in data
  size_t size;
  /// The number of allocated items in data
  size_t capacity;
} catarob_interfaces__srv__ToggleProps_Response__Sequence;

#ifdef __cplusplus
}
#endif

#endif  // CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__STRUCT_H_

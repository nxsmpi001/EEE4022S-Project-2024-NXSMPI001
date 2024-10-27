// generated from rosidl_typesupport_introspection_c/resource/idl__type_support.c.em
// with input from catarob_interfaces:msg\StateEstimate.idl
// generated code does not contain a copyright notice

#include <stddef.h>
#include "catarob_interfaces/msg/detail/state_estimate__rosidl_typesupport_introspection_c.h"
#include "catarob_interfaces/msg/rosidl_typesupport_introspection_c__visibility_control.h"
#include "rosidl_typesupport_introspection_c/field_types.h"
#include "rosidl_typesupport_introspection_c/identifier.h"
#include "rosidl_typesupport_introspection_c/message_introspection.h"
#include "catarob_interfaces/msg/detail/state_estimate__functions.h"
#include "catarob_interfaces/msg/detail/state_estimate__struct.h"


#ifdef __cplusplus
extern "C"
{
#endif

void catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_init_function(
  void * message_memory, enum rosidl_runtime_c__message_initialization _init)
{
  // TODO(karsten1987): initializers are not yet implemented for typesupport c
  // see https://github.com/ros2/ros2/issues/397
  (void) _init;
  catarob_interfaces__msg__StateEstimate__init(message_memory);
}

void catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_fini_function(void * message_memory)
{
  catarob_interfaces__msg__StateEstimate__fini(message_memory);
}

size_t catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__size_function__StateEstimate__covariance(
  const void * untyped_member)
{
  (void)untyped_member;
  return 36;
}

const void * catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__get_const_function__StateEstimate__covariance(
  const void * untyped_member, size_t index)
{
  const double * member =
    (const double *)(untyped_member);
  return &member[index];
}

void * catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__get_function__StateEstimate__covariance(
  void * untyped_member, size_t index)
{
  double * member =
    (double *)(untyped_member);
  return &member[index];
}

void catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__fetch_function__StateEstimate__covariance(
  const void * untyped_member, size_t index, void * untyped_value)
{
  const double * item =
    ((const double *)
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__get_const_function__StateEstimate__covariance(untyped_member, index));
  double * value =
    (double *)(untyped_value);
  *value = *item;
}

void catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__assign_function__StateEstimate__covariance(
  void * untyped_member, size_t index, const void * untyped_value)
{
  double * item =
    ((double *)
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__get_function__StateEstimate__covariance(untyped_member, index));
  const double * value =
    (const double *)(untyped_value);
  *item = *value;
}

static rosidl_typesupport_introspection_c__MessageMember catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_member_array[7] = {
  {
    "x",  // name
    rosidl_typesupport_introspection_c__ROS_TYPE_DOUBLE,  // type
    0,  // upper bound of string
    NULL,  // members of sub message
    false,  // is array
    0,  // array size
    false,  // is upper bound
    offsetof(catarob_interfaces__msg__StateEstimate, x),  // bytes offset in struct
    NULL,  // default value
    NULL,  // size() function pointer
    NULL,  // get_const(index) function pointer
    NULL,  // get(index) function pointer
    NULL,  // fetch(index, &value) function pointer
    NULL,  // assign(index, value) function pointer
    NULL  // resize(index) function pointer
  },
  {
    "y",  // name
    rosidl_typesupport_introspection_c__ROS_TYPE_DOUBLE,  // type
    0,  // upper bound of string
    NULL,  // members of sub message
    false,  // is array
    0,  // array size
    false,  // is upper bound
    offsetof(catarob_interfaces__msg__StateEstimate, y),  // bytes offset in struct
    NULL,  // default value
    NULL,  // size() function pointer
    NULL,  // get_const(index) function pointer
    NULL,  // get(index) function pointer
    NULL,  // fetch(index, &value) function pointer
    NULL,  // assign(index, value) function pointer
    NULL  // resize(index) function pointer
  },
  {
    "psi",  // name
    rosidl_typesupport_introspection_c__ROS_TYPE_DOUBLE,  // type
    0,  // upper bound of string
    NULL,  // members of sub message
    false,  // is array
    0,  // array size
    false,  // is upper bound
    offsetof(catarob_interfaces__msg__StateEstimate, psi),  // bytes offset in struct
    NULL,  // default value
    NULL,  // size() function pointer
    NULL,  // get_const(index) function pointer
    NULL,  // get(index) function pointer
    NULL,  // fetch(index, &value) function pointer
    NULL,  // assign(index, value) function pointer
    NULL  // resize(index) function pointer
  },
  {
    "u",  // name
    rosidl_typesupport_introspection_c__ROS_TYPE_DOUBLE,  // type
    0,  // upper bound of string
    NULL,  // members of sub message
    false,  // is array
    0,  // array size
    false,  // is upper bound
    offsetof(catarob_interfaces__msg__StateEstimate, u),  // bytes offset in struct
    NULL,  // default value
    NULL,  // size() function pointer
    NULL,  // get_const(index) function pointer
    NULL,  // get(index) function pointer
    NULL,  // fetch(index, &value) function pointer
    NULL,  // assign(index, value) function pointer
    NULL  // resize(index) function pointer
  },
  {
    "v",  // name
    rosidl_typesupport_introspection_c__ROS_TYPE_DOUBLE,  // type
    0,  // upper bound of string
    NULL,  // members of sub message
    false,  // is array
    0,  // array size
    false,  // is upper bound
    offsetof(catarob_interfaces__msg__StateEstimate, v),  // bytes offset in struct
    NULL,  // default value
    NULL,  // size() function pointer
    NULL,  // get_const(index) function pointer
    NULL,  // get(index) function pointer
    NULL,  // fetch(index, &value) function pointer
    NULL,  // assign(index, value) function pointer
    NULL  // resize(index) function pointer
  },
  {
    "r",  // name
    rosidl_typesupport_introspection_c__ROS_TYPE_DOUBLE,  // type
    0,  // upper bound of string
    NULL,  // members of sub message
    false,  // is array
    0,  // array size
    false,  // is upper bound
    offsetof(catarob_interfaces__msg__StateEstimate, r),  // bytes offset in struct
    NULL,  // default value
    NULL,  // size() function pointer
    NULL,  // get_const(index) function pointer
    NULL,  // get(index) function pointer
    NULL,  // fetch(index, &value) function pointer
    NULL,  // assign(index, value) function pointer
    NULL  // resize(index) function pointer
  },
  {
    "covariance",  // name
    rosidl_typesupport_introspection_c__ROS_TYPE_DOUBLE,  // type
    0,  // upper bound of string
    NULL,  // members of sub message
    true,  // is array
    36,  // array size
    false,  // is upper bound
    offsetof(catarob_interfaces__msg__StateEstimate, covariance),  // bytes offset in struct
    NULL,  // default value
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__size_function__StateEstimate__covariance,  // size() function pointer
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__get_const_function__StateEstimate__covariance,  // get_const(index) function pointer
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__get_function__StateEstimate__covariance,  // get(index) function pointer
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__fetch_function__StateEstimate__covariance,  // fetch(index, &value) function pointer
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__assign_function__StateEstimate__covariance,  // assign(index, value) function pointer
    NULL  // resize(index) function pointer
  }
};

static const rosidl_typesupport_introspection_c__MessageMembers catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_members = {
  "catarob_interfaces__msg",  // message namespace
  "StateEstimate",  // message name
  7,  // number of fields
  sizeof(catarob_interfaces__msg__StateEstimate),
  catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_member_array,  // message members
  catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_init_function,  // function to initialize message memory (memory has to be allocated)
  catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_fini_function  // function to terminate message instance (will not free memory)
};

// this is not const since it must be initialized on first access
// since C does not allow non-integral compile-time constants
static rosidl_message_type_support_t catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_type_support_handle = {
  0,
  &catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_members,
  get_message_typesupport_handle_function,
};

ROSIDL_TYPESUPPORT_INTROSPECTION_C_EXPORT_catarob_interfaces
const rosidl_message_type_support_t *
ROSIDL_TYPESUPPORT_INTERFACE__MESSAGE_SYMBOL_NAME(rosidl_typesupport_introspection_c, catarob_interfaces, msg, StateEstimate)() {
  if (!catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_type_support_handle.typesupport_identifier) {
    catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_type_support_handle.typesupport_identifier =
      rosidl_typesupport_introspection_c__identifier;
  }
  return &catarob_interfaces__msg__StateEstimate__rosidl_typesupport_introspection_c__StateEstimate_message_type_support_handle;
}
#ifdef __cplusplus
}
#endif

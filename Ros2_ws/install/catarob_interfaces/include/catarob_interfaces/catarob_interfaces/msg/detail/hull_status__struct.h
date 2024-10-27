// generated from rosidl_generator_c/resource/idl__struct.h.em
// with input from catarob_interfaces:msg/HullStatus.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__STRUCT_H_
#define CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__STRUCT_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>


// Constants defined in the message

// Include directives for member types
// Member 'header'
#include "std_msgs/msg/detail/header__struct.h"

/// Struct defined in msg/HullStatus in the package catarob_interfaces.
typedef struct catarob_interfaces__msg__HullStatus
{
  std_msgs__msg__Header header;
  int32_t status;
  float temp;
  int32_t actual_pwm;
  int32_t actual_pwm_rc;
  int32_t pwm_source;
  float ibat;
  float imot;
  float vbat;
  float actual_imax;
  int32_t tor_rc;
  int32_t water_ingress;
  int32_t pwr_relay;
} catarob_interfaces__msg__HullStatus;

// Struct for a sequence of catarob_interfaces__msg__HullStatus.
typedef struct catarob_interfaces__msg__HullStatus__Sequence
{
  catarob_interfaces__msg__HullStatus * data;
  /// The number of valid items in data
  size_t size;
  /// The number of allocated items in data
  size_t capacity;
} catarob_interfaces__msg__HullStatus__Sequence;

#ifdef __cplusplus
}
#endif

#endif  // CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__STRUCT_H_

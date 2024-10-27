// generated from rosidl_generator_c/resource/idl__struct.h.em
// with input from catarob_interfaces:msg/StateEstimate.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__STRUCT_H_
#define CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__STRUCT_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>


// Constants defined in the message

/// Struct defined in msg/StateEstimate in the package catarob_interfaces.
typedef struct catarob_interfaces__msg__StateEstimate
{
  double x;
  double y;
  double psi;
  double u;
  double v;
  double r;
  double covariance[36];
} catarob_interfaces__msg__StateEstimate;

// Struct for a sequence of catarob_interfaces__msg__StateEstimate.
typedef struct catarob_interfaces__msg__StateEstimate__Sequence
{
  catarob_interfaces__msg__StateEstimate * data;
  /// The number of valid items in data
  size_t size;
  /// The number of allocated items in data
  size_t capacity;
} catarob_interfaces__msg__StateEstimate__Sequence;

#ifdef __cplusplus
}
#endif

#endif  // CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__STRUCT_H_

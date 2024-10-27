// generated from rosidl_generator_c/resource/idl__functions.h.em
// with input from catarob_interfaces:msg/StateEstimate.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__FUNCTIONS_H_
#define CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__FUNCTIONS_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>
#include <stdlib.h>

#include "rosidl_runtime_c/visibility_control.h"
#include "catarob_interfaces/msg/rosidl_generator_c__visibility_control.h"

#include "catarob_interfaces/msg/detail/state_estimate__struct.h"

/// Initialize msg/StateEstimate message.
/**
 * If the init function is called twice for the same message without
 * calling fini inbetween previously allocated memory will be leaked.
 * \param[in,out] msg The previously allocated message pointer.
 * Fields without a default value will not be initialized by this function.
 * You might want to call memset(msg, 0, sizeof(
 * catarob_interfaces__msg__StateEstimate
 * )) before or use
 * catarob_interfaces__msg__StateEstimate__create()
 * to allocate and initialize the message.
 * \return true if initialization was successful, otherwise false
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
bool
catarob_interfaces__msg__StateEstimate__init(catarob_interfaces__msg__StateEstimate * msg);

/// Finalize msg/StateEstimate message.
/**
 * \param[in,out] msg The allocated message pointer.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
void
catarob_interfaces__msg__StateEstimate__fini(catarob_interfaces__msg__StateEstimate * msg);

/// Create msg/StateEstimate message.
/**
 * It allocates the memory for the message, sets the memory to zero, and
 * calls
 * catarob_interfaces__msg__StateEstimate__init().
 * \return The pointer to the initialized message if successful,
 * otherwise NULL
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
catarob_interfaces__msg__StateEstimate *
catarob_interfaces__msg__StateEstimate__create();

/// Destroy msg/StateEstimate message.
/**
 * It calls
 * catarob_interfaces__msg__StateEstimate__fini()
 * and frees the memory of the message.
 * \param[in,out] msg The allocated message pointer.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
void
catarob_interfaces__msg__StateEstimate__destroy(catarob_interfaces__msg__StateEstimate * msg);

/// Check for msg/StateEstimate message equality.
/**
 * \param[in] lhs The message on the left hand size of the equality operator.
 * \param[in] rhs The message on the right hand size of the equality operator.
 * \return true if messages are equal, otherwise false.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
bool
catarob_interfaces__msg__StateEstimate__are_equal(const catarob_interfaces__msg__StateEstimate * lhs, const catarob_interfaces__msg__StateEstimate * rhs);

/// Copy a msg/StateEstimate message.
/**
 * This functions performs a deep copy, as opposed to the shallow copy that
 * plain assignment yields.
 *
 * \param[in] input The source message pointer.
 * \param[out] output The target message pointer, which must
 *   have been initialized before calling this function.
 * \return true if successful, or false if either pointer is null
 *   or memory allocation fails.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
bool
catarob_interfaces__msg__StateEstimate__copy(
  const catarob_interfaces__msg__StateEstimate * input,
  catarob_interfaces__msg__StateEstimate * output);

/// Initialize array of msg/StateEstimate messages.
/**
 * It allocates the memory for the number of elements and calls
 * catarob_interfaces__msg__StateEstimate__init()
 * for each element of the array.
 * \param[in,out] array The allocated array pointer.
 * \param[in] size The size / capacity of the array.
 * \return true if initialization was successful, otherwise false
 * If the array pointer is valid and the size is zero it is guaranteed
 # to return true.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
bool
catarob_interfaces__msg__StateEstimate__Sequence__init(catarob_interfaces__msg__StateEstimate__Sequence * array, size_t size);

/// Finalize array of msg/StateEstimate messages.
/**
 * It calls
 * catarob_interfaces__msg__StateEstimate__fini()
 * for each element of the array and frees the memory for the number of
 * elements.
 * \param[in,out] array The initialized array pointer.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
void
catarob_interfaces__msg__StateEstimate__Sequence__fini(catarob_interfaces__msg__StateEstimate__Sequence * array);

/// Create array of msg/StateEstimate messages.
/**
 * It allocates the memory for the array and calls
 * catarob_interfaces__msg__StateEstimate__Sequence__init().
 * \param[in] size The size / capacity of the array.
 * \return The pointer to the initialized array if successful, otherwise NULL
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
catarob_interfaces__msg__StateEstimate__Sequence *
catarob_interfaces__msg__StateEstimate__Sequence__create(size_t size);

/// Destroy array of msg/StateEstimate messages.
/**
 * It calls
 * catarob_interfaces__msg__StateEstimate__Sequence__fini()
 * on the array,
 * and frees the memory of the array.
 * \param[in,out] array The initialized array pointer.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
void
catarob_interfaces__msg__StateEstimate__Sequence__destroy(catarob_interfaces__msg__StateEstimate__Sequence * array);

/// Check for msg/StateEstimate message array equality.
/**
 * \param[in] lhs The message array on the left hand size of the equality operator.
 * \param[in] rhs The message array on the right hand size of the equality operator.
 * \return true if message arrays are equal in size and content, otherwise false.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
bool
catarob_interfaces__msg__StateEstimate__Sequence__are_equal(const catarob_interfaces__msg__StateEstimate__Sequence * lhs, const catarob_interfaces__msg__StateEstimate__Sequence * rhs);

/// Copy an array of msg/StateEstimate messages.
/**
 * This functions performs a deep copy, as opposed to the shallow copy that
 * plain assignment yields.
 *
 * \param[in] input The source array pointer.
 * \param[out] output The target array pointer, which must
 *   have been initialized before calling this function.
 * \return true if successful, or false if either pointer
 *   is null or memory allocation fails.
 */
ROSIDL_GENERATOR_C_PUBLIC_catarob_interfaces
bool
catarob_interfaces__msg__StateEstimate__Sequence__copy(
  const catarob_interfaces__msg__StateEstimate__Sequence * input,
  catarob_interfaces__msg__StateEstimate__Sequence * output);

#ifdef __cplusplus
}
#endif

#endif  // CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__FUNCTIONS_H_

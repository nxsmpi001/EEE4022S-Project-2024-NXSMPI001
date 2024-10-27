// generated from rosidl_generator_c/resource/idl__functions.c.em
// with input from catarob_interfaces:msg/StateEstimate.idl
// generated code does not contain a copyright notice
#include "catarob_interfaces/msg/detail/state_estimate__functions.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "rcutils/allocator.h"


bool
catarob_interfaces__msg__StateEstimate__init(catarob_interfaces__msg__StateEstimate * msg)
{
  if (!msg) {
    return false;
  }
  // x
  // y
  // psi
  // u
  // v
  // r
  // covariance
  return true;
}

void
catarob_interfaces__msg__StateEstimate__fini(catarob_interfaces__msg__StateEstimate * msg)
{
  if (!msg) {
    return;
  }
  // x
  // y
  // psi
  // u
  // v
  // r
  // covariance
}

bool
catarob_interfaces__msg__StateEstimate__are_equal(const catarob_interfaces__msg__StateEstimate * lhs, const catarob_interfaces__msg__StateEstimate * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  // x
  if (lhs->x != rhs->x) {
    return false;
  }
  // y
  if (lhs->y != rhs->y) {
    return false;
  }
  // psi
  if (lhs->psi != rhs->psi) {
    return false;
  }
  // u
  if (lhs->u != rhs->u) {
    return false;
  }
  // v
  if (lhs->v != rhs->v) {
    return false;
  }
  // r
  if (lhs->r != rhs->r) {
    return false;
  }
  // covariance
  for (size_t i = 0; i < 36; ++i) {
    if (lhs->covariance[i] != rhs->covariance[i]) {
      return false;
    }
  }
  return true;
}

bool
catarob_interfaces__msg__StateEstimate__copy(
  const catarob_interfaces__msg__StateEstimate * input,
  catarob_interfaces__msg__StateEstimate * output)
{
  if (!input || !output) {
    return false;
  }
  // x
  output->x = input->x;
  // y
  output->y = input->y;
  // psi
  output->psi = input->psi;
  // u
  output->u = input->u;
  // v
  output->v = input->v;
  // r
  output->r = input->r;
  // covariance
  for (size_t i = 0; i < 36; ++i) {
    output->covariance[i] = input->covariance[i];
  }
  return true;
}

catarob_interfaces__msg__StateEstimate *
catarob_interfaces__msg__StateEstimate__create()
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__msg__StateEstimate * msg = (catarob_interfaces__msg__StateEstimate *)allocator.allocate(sizeof(catarob_interfaces__msg__StateEstimate), allocator.state);
  if (!msg) {
    return NULL;
  }
  memset(msg, 0, sizeof(catarob_interfaces__msg__StateEstimate));
  bool success = catarob_interfaces__msg__StateEstimate__init(msg);
  if (!success) {
    allocator.deallocate(msg, allocator.state);
    return NULL;
  }
  return msg;
}

void
catarob_interfaces__msg__StateEstimate__destroy(catarob_interfaces__msg__StateEstimate * msg)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (msg) {
    catarob_interfaces__msg__StateEstimate__fini(msg);
  }
  allocator.deallocate(msg, allocator.state);
}


bool
catarob_interfaces__msg__StateEstimate__Sequence__init(catarob_interfaces__msg__StateEstimate__Sequence * array, size_t size)
{
  if (!array) {
    return false;
  }
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__msg__StateEstimate * data = NULL;

  if (size) {
    data = (catarob_interfaces__msg__StateEstimate *)allocator.zero_allocate(size, sizeof(catarob_interfaces__msg__StateEstimate), allocator.state);
    if (!data) {
      return false;
    }
    // initialize all array elements
    size_t i;
    for (i = 0; i < size; ++i) {
      bool success = catarob_interfaces__msg__StateEstimate__init(&data[i]);
      if (!success) {
        break;
      }
    }
    if (i < size) {
      // if initialization failed finalize the already initialized array elements
      for (; i > 0; --i) {
        catarob_interfaces__msg__StateEstimate__fini(&data[i - 1]);
      }
      allocator.deallocate(data, allocator.state);
      return false;
    }
  }
  array->data = data;
  array->size = size;
  array->capacity = size;
  return true;
}

void
catarob_interfaces__msg__StateEstimate__Sequence__fini(catarob_interfaces__msg__StateEstimate__Sequence * array)
{
  if (!array) {
    return;
  }
  rcutils_allocator_t allocator = rcutils_get_default_allocator();

  if (array->data) {
    // ensure that data and capacity values are consistent
    assert(array->capacity > 0);
    // finalize all array elements
    for (size_t i = 0; i < array->capacity; ++i) {
      catarob_interfaces__msg__StateEstimate__fini(&array->data[i]);
    }
    allocator.deallocate(array->data, allocator.state);
    array->data = NULL;
    array->size = 0;
    array->capacity = 0;
  } else {
    // ensure that data, size, and capacity values are consistent
    assert(0 == array->size);
    assert(0 == array->capacity);
  }
}

catarob_interfaces__msg__StateEstimate__Sequence *
catarob_interfaces__msg__StateEstimate__Sequence__create(size_t size)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__msg__StateEstimate__Sequence * array = (catarob_interfaces__msg__StateEstimate__Sequence *)allocator.allocate(sizeof(catarob_interfaces__msg__StateEstimate__Sequence), allocator.state);
  if (!array) {
    return NULL;
  }
  bool success = catarob_interfaces__msg__StateEstimate__Sequence__init(array, size);
  if (!success) {
    allocator.deallocate(array, allocator.state);
    return NULL;
  }
  return array;
}

void
catarob_interfaces__msg__StateEstimate__Sequence__destroy(catarob_interfaces__msg__StateEstimate__Sequence * array)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (array) {
    catarob_interfaces__msg__StateEstimate__Sequence__fini(array);
  }
  allocator.deallocate(array, allocator.state);
}

bool
catarob_interfaces__msg__StateEstimate__Sequence__are_equal(const catarob_interfaces__msg__StateEstimate__Sequence * lhs, const catarob_interfaces__msg__StateEstimate__Sequence * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  if (lhs->size != rhs->size) {
    return false;
  }
  for (size_t i = 0; i < lhs->size; ++i) {
    if (!catarob_interfaces__msg__StateEstimate__are_equal(&(lhs->data[i]), &(rhs->data[i]))) {
      return false;
    }
  }
  return true;
}

bool
catarob_interfaces__msg__StateEstimate__Sequence__copy(
  const catarob_interfaces__msg__StateEstimate__Sequence * input,
  catarob_interfaces__msg__StateEstimate__Sequence * output)
{
  if (!input || !output) {
    return false;
  }
  if (output->capacity < input->size) {
    const size_t allocation_size =
      input->size * sizeof(catarob_interfaces__msg__StateEstimate);
    rcutils_allocator_t allocator = rcutils_get_default_allocator();
    catarob_interfaces__msg__StateEstimate * data =
      (catarob_interfaces__msg__StateEstimate *)allocator.reallocate(
      output->data, allocation_size, allocator.state);
    if (!data) {
      return false;
    }
    // If reallocation succeeded, memory may or may not have been moved
    // to fulfill the allocation request, invalidating output->data.
    output->data = data;
    for (size_t i = output->capacity; i < input->size; ++i) {
      if (!catarob_interfaces__msg__StateEstimate__init(&output->data[i])) {
        // If initialization of any new item fails, roll back
        // all previously initialized items. Existing items
        // in output are to be left unmodified.
        for (; i-- > output->capacity; ) {
          catarob_interfaces__msg__StateEstimate__fini(&output->data[i]);
        }
        return false;
      }
    }
    output->capacity = input->size;
  }
  output->size = input->size;
  for (size_t i = 0; i < input->size; ++i) {
    if (!catarob_interfaces__msg__StateEstimate__copy(
        &(input->data[i]), &(output->data[i])))
    {
      return false;
    }
  }
  return true;
}

// generated from rosidl_generator_c/resource/idl__functions.c.em
// with input from catarob_interfaces:msg/HullStatus.idl
// generated code does not contain a copyright notice
#include "catarob_interfaces/msg/detail/hull_status__functions.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "rcutils/allocator.h"


// Include directives for member types
// Member `header`
#include "std_msgs/msg/detail/header__functions.h"

bool
catarob_interfaces__msg__HullStatus__init(catarob_interfaces__msg__HullStatus * msg)
{
  if (!msg) {
    return false;
  }
  // header
  if (!std_msgs__msg__Header__init(&msg->header)) {
    catarob_interfaces__msg__HullStatus__fini(msg);
    return false;
  }
  // status
  // temp
  // actual_pwm
  // actual_pwm_rc
  // pwm_source
  // ibat
  // imot
  // vbat
  // actual_imax
  // tor_rc
  // water_ingress
  // pwr_relay
  return true;
}

void
catarob_interfaces__msg__HullStatus__fini(catarob_interfaces__msg__HullStatus * msg)
{
  if (!msg) {
    return;
  }
  // header
  std_msgs__msg__Header__fini(&msg->header);
  // status
  // temp
  // actual_pwm
  // actual_pwm_rc
  // pwm_source
  // ibat
  // imot
  // vbat
  // actual_imax
  // tor_rc
  // water_ingress
  // pwr_relay
}

bool
catarob_interfaces__msg__HullStatus__are_equal(const catarob_interfaces__msg__HullStatus * lhs, const catarob_interfaces__msg__HullStatus * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  // header
  if (!std_msgs__msg__Header__are_equal(
      &(lhs->header), &(rhs->header)))
  {
    return false;
  }
  // status
  if (lhs->status != rhs->status) {
    return false;
  }
  // temp
  if (lhs->temp != rhs->temp) {
    return false;
  }
  // actual_pwm
  if (lhs->actual_pwm != rhs->actual_pwm) {
    return false;
  }
  // actual_pwm_rc
  if (lhs->actual_pwm_rc != rhs->actual_pwm_rc) {
    return false;
  }
  // pwm_source
  if (lhs->pwm_source != rhs->pwm_source) {
    return false;
  }
  // ibat
  if (lhs->ibat != rhs->ibat) {
    return false;
  }
  // imot
  if (lhs->imot != rhs->imot) {
    return false;
  }
  // vbat
  if (lhs->vbat != rhs->vbat) {
    return false;
  }
  // actual_imax
  if (lhs->actual_imax != rhs->actual_imax) {
    return false;
  }
  // tor_rc
  if (lhs->tor_rc != rhs->tor_rc) {
    return false;
  }
  // water_ingress
  if (lhs->water_ingress != rhs->water_ingress) {
    return false;
  }
  // pwr_relay
  if (lhs->pwr_relay != rhs->pwr_relay) {
    return false;
  }
  return true;
}

bool
catarob_interfaces__msg__HullStatus__copy(
  const catarob_interfaces__msg__HullStatus * input,
  catarob_interfaces__msg__HullStatus * output)
{
  if (!input || !output) {
    return false;
  }
  // header
  if (!std_msgs__msg__Header__copy(
      &(input->header), &(output->header)))
  {
    return false;
  }
  // status
  output->status = input->status;
  // temp
  output->temp = input->temp;
  // actual_pwm
  output->actual_pwm = input->actual_pwm;
  // actual_pwm_rc
  output->actual_pwm_rc = input->actual_pwm_rc;
  // pwm_source
  output->pwm_source = input->pwm_source;
  // ibat
  output->ibat = input->ibat;
  // imot
  output->imot = input->imot;
  // vbat
  output->vbat = input->vbat;
  // actual_imax
  output->actual_imax = input->actual_imax;
  // tor_rc
  output->tor_rc = input->tor_rc;
  // water_ingress
  output->water_ingress = input->water_ingress;
  // pwr_relay
  output->pwr_relay = input->pwr_relay;
  return true;
}

catarob_interfaces__msg__HullStatus *
catarob_interfaces__msg__HullStatus__create()
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__msg__HullStatus * msg = (catarob_interfaces__msg__HullStatus *)allocator.allocate(sizeof(catarob_interfaces__msg__HullStatus), allocator.state);
  if (!msg) {
    return NULL;
  }
  memset(msg, 0, sizeof(catarob_interfaces__msg__HullStatus));
  bool success = catarob_interfaces__msg__HullStatus__init(msg);
  if (!success) {
    allocator.deallocate(msg, allocator.state);
    return NULL;
  }
  return msg;
}

void
catarob_interfaces__msg__HullStatus__destroy(catarob_interfaces__msg__HullStatus * msg)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (msg) {
    catarob_interfaces__msg__HullStatus__fini(msg);
  }
  allocator.deallocate(msg, allocator.state);
}


bool
catarob_interfaces__msg__HullStatus__Sequence__init(catarob_interfaces__msg__HullStatus__Sequence * array, size_t size)
{
  if (!array) {
    return false;
  }
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__msg__HullStatus * data = NULL;

  if (size) {
    data = (catarob_interfaces__msg__HullStatus *)allocator.zero_allocate(size, sizeof(catarob_interfaces__msg__HullStatus), allocator.state);
    if (!data) {
      return false;
    }
    // initialize all array elements
    size_t i;
    for (i = 0; i < size; ++i) {
      bool success = catarob_interfaces__msg__HullStatus__init(&data[i]);
      if (!success) {
        break;
      }
    }
    if (i < size) {
      // if initialization failed finalize the already initialized array elements
      for (; i > 0; --i) {
        catarob_interfaces__msg__HullStatus__fini(&data[i - 1]);
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
catarob_interfaces__msg__HullStatus__Sequence__fini(catarob_interfaces__msg__HullStatus__Sequence * array)
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
      catarob_interfaces__msg__HullStatus__fini(&array->data[i]);
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

catarob_interfaces__msg__HullStatus__Sequence *
catarob_interfaces__msg__HullStatus__Sequence__create(size_t size)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__msg__HullStatus__Sequence * array = (catarob_interfaces__msg__HullStatus__Sequence *)allocator.allocate(sizeof(catarob_interfaces__msg__HullStatus__Sequence), allocator.state);
  if (!array) {
    return NULL;
  }
  bool success = catarob_interfaces__msg__HullStatus__Sequence__init(array, size);
  if (!success) {
    allocator.deallocate(array, allocator.state);
    return NULL;
  }
  return array;
}

void
catarob_interfaces__msg__HullStatus__Sequence__destroy(catarob_interfaces__msg__HullStatus__Sequence * array)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (array) {
    catarob_interfaces__msg__HullStatus__Sequence__fini(array);
  }
  allocator.deallocate(array, allocator.state);
}

bool
catarob_interfaces__msg__HullStatus__Sequence__are_equal(const catarob_interfaces__msg__HullStatus__Sequence * lhs, const catarob_interfaces__msg__HullStatus__Sequence * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  if (lhs->size != rhs->size) {
    return false;
  }
  for (size_t i = 0; i < lhs->size; ++i) {
    if (!catarob_interfaces__msg__HullStatus__are_equal(&(lhs->data[i]), &(rhs->data[i]))) {
      return false;
    }
  }
  return true;
}

bool
catarob_interfaces__msg__HullStatus__Sequence__copy(
  const catarob_interfaces__msg__HullStatus__Sequence * input,
  catarob_interfaces__msg__HullStatus__Sequence * output)
{
  if (!input || !output) {
    return false;
  }
  if (output->capacity < input->size) {
    const size_t allocation_size =
      input->size * sizeof(catarob_interfaces__msg__HullStatus);
    rcutils_allocator_t allocator = rcutils_get_default_allocator();
    catarob_interfaces__msg__HullStatus * data =
      (catarob_interfaces__msg__HullStatus *)allocator.reallocate(
      output->data, allocation_size, allocator.state);
    if (!data) {
      return false;
    }
    // If reallocation succeeded, memory may or may not have been moved
    // to fulfill the allocation request, invalidating output->data.
    output->data = data;
    for (size_t i = output->capacity; i < input->size; ++i) {
      if (!catarob_interfaces__msg__HullStatus__init(&output->data[i])) {
        // If initialization of any new item fails, roll back
        // all previously initialized items. Existing items
        // in output are to be left unmodified.
        for (; i-- > output->capacity; ) {
          catarob_interfaces__msg__HullStatus__fini(&output->data[i]);
        }
        return false;
      }
    }
    output->capacity = input->size;
  }
  output->size = input->size;
  for (size_t i = 0; i < input->size; ++i) {
    if (!catarob_interfaces__msg__HullStatus__copy(
        &(input->data[i]), &(output->data[i])))
    {
      return false;
    }
  }
  return true;
}

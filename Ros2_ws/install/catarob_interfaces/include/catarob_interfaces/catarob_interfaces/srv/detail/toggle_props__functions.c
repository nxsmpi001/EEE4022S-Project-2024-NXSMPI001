// generated from rosidl_generator_c/resource/idl__functions.c.em
// with input from catarob_interfaces:srv/ToggleProps.idl
// generated code does not contain a copyright notice
#include "catarob_interfaces/srv/detail/toggle_props__functions.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "rcutils/allocator.h"

bool
catarob_interfaces__srv__ToggleProps_Request__init(catarob_interfaces__srv__ToggleProps_Request * msg)
{
  if (!msg) {
    return false;
  }
  // flag
  return true;
}

void
catarob_interfaces__srv__ToggleProps_Request__fini(catarob_interfaces__srv__ToggleProps_Request * msg)
{
  if (!msg) {
    return;
  }
  // flag
}

bool
catarob_interfaces__srv__ToggleProps_Request__are_equal(const catarob_interfaces__srv__ToggleProps_Request * lhs, const catarob_interfaces__srv__ToggleProps_Request * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  // flag
  if (lhs->flag != rhs->flag) {
    return false;
  }
  return true;
}

bool
catarob_interfaces__srv__ToggleProps_Request__copy(
  const catarob_interfaces__srv__ToggleProps_Request * input,
  catarob_interfaces__srv__ToggleProps_Request * output)
{
  if (!input || !output) {
    return false;
  }
  // flag
  output->flag = input->flag;
  return true;
}

catarob_interfaces__srv__ToggleProps_Request *
catarob_interfaces__srv__ToggleProps_Request__create()
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__srv__ToggleProps_Request * msg = (catarob_interfaces__srv__ToggleProps_Request *)allocator.allocate(sizeof(catarob_interfaces__srv__ToggleProps_Request), allocator.state);
  if (!msg) {
    return NULL;
  }
  memset(msg, 0, sizeof(catarob_interfaces__srv__ToggleProps_Request));
  bool success = catarob_interfaces__srv__ToggleProps_Request__init(msg);
  if (!success) {
    allocator.deallocate(msg, allocator.state);
    return NULL;
  }
  return msg;
}

void
catarob_interfaces__srv__ToggleProps_Request__destroy(catarob_interfaces__srv__ToggleProps_Request * msg)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (msg) {
    catarob_interfaces__srv__ToggleProps_Request__fini(msg);
  }
  allocator.deallocate(msg, allocator.state);
}


bool
catarob_interfaces__srv__ToggleProps_Request__Sequence__init(catarob_interfaces__srv__ToggleProps_Request__Sequence * array, size_t size)
{
  if (!array) {
    return false;
  }
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__srv__ToggleProps_Request * data = NULL;

  if (size) {
    data = (catarob_interfaces__srv__ToggleProps_Request *)allocator.zero_allocate(size, sizeof(catarob_interfaces__srv__ToggleProps_Request), allocator.state);
    if (!data) {
      return false;
    }
    // initialize all array elements
    size_t i;
    for (i = 0; i < size; ++i) {
      bool success = catarob_interfaces__srv__ToggleProps_Request__init(&data[i]);
      if (!success) {
        break;
      }
    }
    if (i < size) {
      // if initialization failed finalize the already initialized array elements
      for (; i > 0; --i) {
        catarob_interfaces__srv__ToggleProps_Request__fini(&data[i - 1]);
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
catarob_interfaces__srv__ToggleProps_Request__Sequence__fini(catarob_interfaces__srv__ToggleProps_Request__Sequence * array)
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
      catarob_interfaces__srv__ToggleProps_Request__fini(&array->data[i]);
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

catarob_interfaces__srv__ToggleProps_Request__Sequence *
catarob_interfaces__srv__ToggleProps_Request__Sequence__create(size_t size)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__srv__ToggleProps_Request__Sequence * array = (catarob_interfaces__srv__ToggleProps_Request__Sequence *)allocator.allocate(sizeof(catarob_interfaces__srv__ToggleProps_Request__Sequence), allocator.state);
  if (!array) {
    return NULL;
  }
  bool success = catarob_interfaces__srv__ToggleProps_Request__Sequence__init(array, size);
  if (!success) {
    allocator.deallocate(array, allocator.state);
    return NULL;
  }
  return array;
}

void
catarob_interfaces__srv__ToggleProps_Request__Sequence__destroy(catarob_interfaces__srv__ToggleProps_Request__Sequence * array)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (array) {
    catarob_interfaces__srv__ToggleProps_Request__Sequence__fini(array);
  }
  allocator.deallocate(array, allocator.state);
}

bool
catarob_interfaces__srv__ToggleProps_Request__Sequence__are_equal(const catarob_interfaces__srv__ToggleProps_Request__Sequence * lhs, const catarob_interfaces__srv__ToggleProps_Request__Sequence * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  if (lhs->size != rhs->size) {
    return false;
  }
  for (size_t i = 0; i < lhs->size; ++i) {
    if (!catarob_interfaces__srv__ToggleProps_Request__are_equal(&(lhs->data[i]), &(rhs->data[i]))) {
      return false;
    }
  }
  return true;
}

bool
catarob_interfaces__srv__ToggleProps_Request__Sequence__copy(
  const catarob_interfaces__srv__ToggleProps_Request__Sequence * input,
  catarob_interfaces__srv__ToggleProps_Request__Sequence * output)
{
  if (!input || !output) {
    return false;
  }
  if (output->capacity < input->size) {
    const size_t allocation_size =
      input->size * sizeof(catarob_interfaces__srv__ToggleProps_Request);
    rcutils_allocator_t allocator = rcutils_get_default_allocator();
    catarob_interfaces__srv__ToggleProps_Request * data =
      (catarob_interfaces__srv__ToggleProps_Request *)allocator.reallocate(
      output->data, allocation_size, allocator.state);
    if (!data) {
      return false;
    }
    // If reallocation succeeded, memory may or may not have been moved
    // to fulfill the allocation request, invalidating output->data.
    output->data = data;
    for (size_t i = output->capacity; i < input->size; ++i) {
      if (!catarob_interfaces__srv__ToggleProps_Request__init(&output->data[i])) {
        // If initialization of any new item fails, roll back
        // all previously initialized items. Existing items
        // in output are to be left unmodified.
        for (; i-- > output->capacity; ) {
          catarob_interfaces__srv__ToggleProps_Request__fini(&output->data[i]);
        }
        return false;
      }
    }
    output->capacity = input->size;
  }
  output->size = input->size;
  for (size_t i = 0; i < input->size; ++i) {
    if (!catarob_interfaces__srv__ToggleProps_Request__copy(
        &(input->data[i]), &(output->data[i])))
    {
      return false;
    }
  }
  return true;
}


bool
catarob_interfaces__srv__ToggleProps_Response__init(catarob_interfaces__srv__ToggleProps_Response * msg)
{
  if (!msg) {
    return false;
  }
  // state
  return true;
}

void
catarob_interfaces__srv__ToggleProps_Response__fini(catarob_interfaces__srv__ToggleProps_Response * msg)
{
  if (!msg) {
    return;
  }
  // state
}

bool
catarob_interfaces__srv__ToggleProps_Response__are_equal(const catarob_interfaces__srv__ToggleProps_Response * lhs, const catarob_interfaces__srv__ToggleProps_Response * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  // state
  if (lhs->state != rhs->state) {
    return false;
  }
  return true;
}

bool
catarob_interfaces__srv__ToggleProps_Response__copy(
  const catarob_interfaces__srv__ToggleProps_Response * input,
  catarob_interfaces__srv__ToggleProps_Response * output)
{
  if (!input || !output) {
    return false;
  }
  // state
  output->state = input->state;
  return true;
}

catarob_interfaces__srv__ToggleProps_Response *
catarob_interfaces__srv__ToggleProps_Response__create()
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__srv__ToggleProps_Response * msg = (catarob_interfaces__srv__ToggleProps_Response *)allocator.allocate(sizeof(catarob_interfaces__srv__ToggleProps_Response), allocator.state);
  if (!msg) {
    return NULL;
  }
  memset(msg, 0, sizeof(catarob_interfaces__srv__ToggleProps_Response));
  bool success = catarob_interfaces__srv__ToggleProps_Response__init(msg);
  if (!success) {
    allocator.deallocate(msg, allocator.state);
    return NULL;
  }
  return msg;
}

void
catarob_interfaces__srv__ToggleProps_Response__destroy(catarob_interfaces__srv__ToggleProps_Response * msg)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (msg) {
    catarob_interfaces__srv__ToggleProps_Response__fini(msg);
  }
  allocator.deallocate(msg, allocator.state);
}


bool
catarob_interfaces__srv__ToggleProps_Response__Sequence__init(catarob_interfaces__srv__ToggleProps_Response__Sequence * array, size_t size)
{
  if (!array) {
    return false;
  }
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__srv__ToggleProps_Response * data = NULL;

  if (size) {
    data = (catarob_interfaces__srv__ToggleProps_Response *)allocator.zero_allocate(size, sizeof(catarob_interfaces__srv__ToggleProps_Response), allocator.state);
    if (!data) {
      return false;
    }
    // initialize all array elements
    size_t i;
    for (i = 0; i < size; ++i) {
      bool success = catarob_interfaces__srv__ToggleProps_Response__init(&data[i]);
      if (!success) {
        break;
      }
    }
    if (i < size) {
      // if initialization failed finalize the already initialized array elements
      for (; i > 0; --i) {
        catarob_interfaces__srv__ToggleProps_Response__fini(&data[i - 1]);
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
catarob_interfaces__srv__ToggleProps_Response__Sequence__fini(catarob_interfaces__srv__ToggleProps_Response__Sequence * array)
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
      catarob_interfaces__srv__ToggleProps_Response__fini(&array->data[i]);
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

catarob_interfaces__srv__ToggleProps_Response__Sequence *
catarob_interfaces__srv__ToggleProps_Response__Sequence__create(size_t size)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  catarob_interfaces__srv__ToggleProps_Response__Sequence * array = (catarob_interfaces__srv__ToggleProps_Response__Sequence *)allocator.allocate(sizeof(catarob_interfaces__srv__ToggleProps_Response__Sequence), allocator.state);
  if (!array) {
    return NULL;
  }
  bool success = catarob_interfaces__srv__ToggleProps_Response__Sequence__init(array, size);
  if (!success) {
    allocator.deallocate(array, allocator.state);
    return NULL;
  }
  return array;
}

void
catarob_interfaces__srv__ToggleProps_Response__Sequence__destroy(catarob_interfaces__srv__ToggleProps_Response__Sequence * array)
{
  rcutils_allocator_t allocator = rcutils_get_default_allocator();
  if (array) {
    catarob_interfaces__srv__ToggleProps_Response__Sequence__fini(array);
  }
  allocator.deallocate(array, allocator.state);
}

bool
catarob_interfaces__srv__ToggleProps_Response__Sequence__are_equal(const catarob_interfaces__srv__ToggleProps_Response__Sequence * lhs, const catarob_interfaces__srv__ToggleProps_Response__Sequence * rhs)
{
  if (!lhs || !rhs) {
    return false;
  }
  if (lhs->size != rhs->size) {
    return false;
  }
  for (size_t i = 0; i < lhs->size; ++i) {
    if (!catarob_interfaces__srv__ToggleProps_Response__are_equal(&(lhs->data[i]), &(rhs->data[i]))) {
      return false;
    }
  }
  return true;
}

bool
catarob_interfaces__srv__ToggleProps_Response__Sequence__copy(
  const catarob_interfaces__srv__ToggleProps_Response__Sequence * input,
  catarob_interfaces__srv__ToggleProps_Response__Sequence * output)
{
  if (!input || !output) {
    return false;
  }
  if (output->capacity < input->size) {
    const size_t allocation_size =
      input->size * sizeof(catarob_interfaces__srv__ToggleProps_Response);
    rcutils_allocator_t allocator = rcutils_get_default_allocator();
    catarob_interfaces__srv__ToggleProps_Response * data =
      (catarob_interfaces__srv__ToggleProps_Response *)allocator.reallocate(
      output->data, allocation_size, allocator.state);
    if (!data) {
      return false;
    }
    // If reallocation succeeded, memory may or may not have been moved
    // to fulfill the allocation request, invalidating output->data.
    output->data = data;
    for (size_t i = output->capacity; i < input->size; ++i) {
      if (!catarob_interfaces__srv__ToggleProps_Response__init(&output->data[i])) {
        // If initialization of any new item fails, roll back
        // all previously initialized items. Existing items
        // in output are to be left unmodified.
        for (; i-- > output->capacity; ) {
          catarob_interfaces__srv__ToggleProps_Response__fini(&output->data[i]);
        }
        return false;
      }
    }
    output->capacity = input->size;
  }
  output->size = input->size;
  for (size_t i = 0; i < input->size; ++i) {
    if (!catarob_interfaces__srv__ToggleProps_Response__copy(
        &(input->data[i]), &(output->data[i])))
    {
      return false;
    }
  }
  return true;
}

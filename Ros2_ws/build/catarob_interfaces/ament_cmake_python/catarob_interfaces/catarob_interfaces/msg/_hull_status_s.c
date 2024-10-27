// generated from rosidl_generator_py/resource/_idl_support.c.em
// with input from catarob_interfaces:msg/HullStatus.idl
// generated code does not contain a copyright notice
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdbool.h>
#ifndef _WIN32
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wunused-function"
#endif
#include "numpy/ndarrayobject.h"
#ifndef _WIN32
# pragma GCC diagnostic pop
#endif
#include "rosidl_runtime_c/visibility_control.h"
#include "catarob_interfaces/msg/detail/hull_status__struct.h"
#include "catarob_interfaces/msg/detail/hull_status__functions.h"

ROSIDL_GENERATOR_C_IMPORT
bool std_msgs__msg__header__convert_from_py(PyObject * _pymsg, void * _ros_message);
ROSIDL_GENERATOR_C_IMPORT
PyObject * std_msgs__msg__header__convert_to_py(void * raw_ros_message);

ROSIDL_GENERATOR_C_EXPORT
bool catarob_interfaces__msg__hull_status__convert_from_py(PyObject * _pymsg, void * _ros_message)
{
  // check that the passed message is of the expected Python class
  {
    char full_classname_dest[47];
    {
      char * class_name = NULL;
      char * module_name = NULL;
      {
        PyObject * class_attr = PyObject_GetAttrString(_pymsg, "__class__");
        if (class_attr) {
          PyObject * name_attr = PyObject_GetAttrString(class_attr, "__name__");
          if (name_attr) {
            class_name = (char *)PyUnicode_1BYTE_DATA(name_attr);
            Py_DECREF(name_attr);
          }
          PyObject * module_attr = PyObject_GetAttrString(class_attr, "__module__");
          if (module_attr) {
            module_name = (char *)PyUnicode_1BYTE_DATA(module_attr);
            Py_DECREF(module_attr);
          }
          Py_DECREF(class_attr);
        }
      }
      if (!class_name || !module_name) {
        return false;
      }
      snprintf(full_classname_dest, sizeof(full_classname_dest), "%s.%s", module_name, class_name);
    }
    assert(strncmp("catarob_interfaces.msg._hull_status.HullStatus", full_classname_dest, 46) == 0);
  }
  catarob_interfaces__msg__HullStatus * ros_message = _ros_message;
  {  // header
    PyObject * field = PyObject_GetAttrString(_pymsg, "header");
    if (!field) {
      return false;
    }
    if (!std_msgs__msg__header__convert_from_py(field, &ros_message->header)) {
      Py_DECREF(field);
      return false;
    }
    Py_DECREF(field);
  }
  {  // status
    PyObject * field = PyObject_GetAttrString(_pymsg, "status");
    if (!field) {
      return false;
    }
    assert(PyLong_Check(field));
    ros_message->status = (int32_t)PyLong_AsLong(field);
    Py_DECREF(field);
  }
  {  // temp
    PyObject * field = PyObject_GetAttrString(_pymsg, "temp");
    if (!field) {
      return false;
    }
    assert(PyFloat_Check(field));
    ros_message->temp = (float)PyFloat_AS_DOUBLE(field);
    Py_DECREF(field);
  }
  {  // actual_pwm
    PyObject * field = PyObject_GetAttrString(_pymsg, "actual_pwm");
    if (!field) {
      return false;
    }
    assert(PyLong_Check(field));
    ros_message->actual_pwm = (int32_t)PyLong_AsLong(field);
    Py_DECREF(field);
  }
  {  // actual_pwm_rc
    PyObject * field = PyObject_GetAttrString(_pymsg, "actual_pwm_rc");
    if (!field) {
      return false;
    }
    assert(PyLong_Check(field));
    ros_message->actual_pwm_rc = (int32_t)PyLong_AsLong(field);
    Py_DECREF(field);
  }
  {  // pwm_source
    PyObject * field = PyObject_GetAttrString(_pymsg, "pwm_source");
    if (!field) {
      return false;
    }
    assert(PyLong_Check(field));
    ros_message->pwm_source = (int32_t)PyLong_AsLong(field);
    Py_DECREF(field);
  }
  {  // ibat
    PyObject * field = PyObject_GetAttrString(_pymsg, "ibat");
    if (!field) {
      return false;
    }
    assert(PyFloat_Check(field));
    ros_message->ibat = (float)PyFloat_AS_DOUBLE(field);
    Py_DECREF(field);
  }
  {  // imot
    PyObject * field = PyObject_GetAttrString(_pymsg, "imot");
    if (!field) {
      return false;
    }
    assert(PyFloat_Check(field));
    ros_message->imot = (float)PyFloat_AS_DOUBLE(field);
    Py_DECREF(field);
  }
  {  // vbat
    PyObject * field = PyObject_GetAttrString(_pymsg, "vbat");
    if (!field) {
      return false;
    }
    assert(PyFloat_Check(field));
    ros_message->vbat = (float)PyFloat_AS_DOUBLE(field);
    Py_DECREF(field);
  }
  {  // actual_imax
    PyObject * field = PyObject_GetAttrString(_pymsg, "actual_imax");
    if (!field) {
      return false;
    }
    assert(PyFloat_Check(field));
    ros_message->actual_imax = (float)PyFloat_AS_DOUBLE(field);
    Py_DECREF(field);
  }
  {  // tor_rc
    PyObject * field = PyObject_GetAttrString(_pymsg, "tor_rc");
    if (!field) {
      return false;
    }
    assert(PyLong_Check(field));
    ros_message->tor_rc = (int32_t)PyLong_AsLong(field);
    Py_DECREF(field);
  }
  {  // water_ingress
    PyObject * field = PyObject_GetAttrString(_pymsg, "water_ingress");
    if (!field) {
      return false;
    }
    assert(PyLong_Check(field));
    ros_message->water_ingress = (int32_t)PyLong_AsLong(field);
    Py_DECREF(field);
  }
  {  // pwr_relay
    PyObject * field = PyObject_GetAttrString(_pymsg, "pwr_relay");
    if (!field) {
      return false;
    }
    assert(PyLong_Check(field));
    ros_message->pwr_relay = (int32_t)PyLong_AsLong(field);
    Py_DECREF(field);
  }

  return true;
}

ROSIDL_GENERATOR_C_EXPORT
PyObject * catarob_interfaces__msg__hull_status__convert_to_py(void * raw_ros_message)
{
  /* NOTE(esteve): Call constructor of HullStatus */
  PyObject * _pymessage = NULL;
  {
    PyObject * pymessage_module = PyImport_ImportModule("catarob_interfaces.msg._hull_status");
    assert(pymessage_module);
    PyObject * pymessage_class = PyObject_GetAttrString(pymessage_module, "HullStatus");
    assert(pymessage_class);
    Py_DECREF(pymessage_module);
    _pymessage = PyObject_CallObject(pymessage_class, NULL);
    Py_DECREF(pymessage_class);
    if (!_pymessage) {
      return NULL;
    }
  }
  catarob_interfaces__msg__HullStatus * ros_message = (catarob_interfaces__msg__HullStatus *)raw_ros_message;
  {  // header
    PyObject * field = NULL;
    field = std_msgs__msg__header__convert_to_py(&ros_message->header);
    if (!field) {
      return NULL;
    }
    {
      int rc = PyObject_SetAttrString(_pymessage, "header", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // status
    PyObject * field = NULL;
    field = PyLong_FromLong(ros_message->status);
    {
      int rc = PyObject_SetAttrString(_pymessage, "status", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // temp
    PyObject * field = NULL;
    field = PyFloat_FromDouble(ros_message->temp);
    {
      int rc = PyObject_SetAttrString(_pymessage, "temp", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // actual_pwm
    PyObject * field = NULL;
    field = PyLong_FromLong(ros_message->actual_pwm);
    {
      int rc = PyObject_SetAttrString(_pymessage, "actual_pwm", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // actual_pwm_rc
    PyObject * field = NULL;
    field = PyLong_FromLong(ros_message->actual_pwm_rc);
    {
      int rc = PyObject_SetAttrString(_pymessage, "actual_pwm_rc", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // pwm_source
    PyObject * field = NULL;
    field = PyLong_FromLong(ros_message->pwm_source);
    {
      int rc = PyObject_SetAttrString(_pymessage, "pwm_source", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // ibat
    PyObject * field = NULL;
    field = PyFloat_FromDouble(ros_message->ibat);
    {
      int rc = PyObject_SetAttrString(_pymessage, "ibat", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // imot
    PyObject * field = NULL;
    field = PyFloat_FromDouble(ros_message->imot);
    {
      int rc = PyObject_SetAttrString(_pymessage, "imot", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // vbat
    PyObject * field = NULL;
    field = PyFloat_FromDouble(ros_message->vbat);
    {
      int rc = PyObject_SetAttrString(_pymessage, "vbat", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // actual_imax
    PyObject * field = NULL;
    field = PyFloat_FromDouble(ros_message->actual_imax);
    {
      int rc = PyObject_SetAttrString(_pymessage, "actual_imax", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // tor_rc
    PyObject * field = NULL;
    field = PyLong_FromLong(ros_message->tor_rc);
    {
      int rc = PyObject_SetAttrString(_pymessage, "tor_rc", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // water_ingress
    PyObject * field = NULL;
    field = PyLong_FromLong(ros_message->water_ingress);
    {
      int rc = PyObject_SetAttrString(_pymessage, "water_ingress", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }
  {  // pwr_relay
    PyObject * field = NULL;
    field = PyLong_FromLong(ros_message->pwr_relay);
    {
      int rc = PyObject_SetAttrString(_pymessage, "pwr_relay", field);
      Py_DECREF(field);
      if (rc) {
        return NULL;
      }
    }
  }

  // ownership of _pymessage is transferred to the caller
  return _pymessage;
}

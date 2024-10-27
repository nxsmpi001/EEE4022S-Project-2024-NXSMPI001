// generated from rosidl_generator_cpp/resource/idl__struct.hpp.em
// with input from catarob_interfaces:msg\HullStatus.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__STRUCT_HPP_
#define CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__STRUCT_HPP_

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "rosidl_runtime_cpp/bounded_vector.hpp"
#include "rosidl_runtime_cpp/message_initialization.hpp"


// Include directives for member types
// Member 'header'
#include "std_msgs/msg/detail/header__struct.hpp"

#ifndef _WIN32
# define DEPRECATED__catarob_interfaces__msg__HullStatus __attribute__((deprecated))
#else
# define DEPRECATED__catarob_interfaces__msg__HullStatus __declspec(deprecated)
#endif

namespace catarob_interfaces
{

namespace msg
{

// message struct
template<class ContainerAllocator>
struct HullStatus_
{
  using Type = HullStatus_<ContainerAllocator>;

  explicit HullStatus_(rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  : header(_init)
  {
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->status = 0l;
      this->temp = 0.0f;
      this->actual_pwm = 0l;
      this->actual_pwm_rc = 0l;
      this->pwm_source = 0l;
      this->ibat = 0.0f;
      this->imot = 0.0f;
      this->vbat = 0.0f;
      this->actual_imax = 0.0f;
      this->tor_rc = 0l;
      this->water_ingress = 0l;
      this->pwr_relay = 0l;
    }
  }

  explicit HullStatus_(const ContainerAllocator & _alloc, rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  : header(_alloc, _init)
  {
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->status = 0l;
      this->temp = 0.0f;
      this->actual_pwm = 0l;
      this->actual_pwm_rc = 0l;
      this->pwm_source = 0l;
      this->ibat = 0.0f;
      this->imot = 0.0f;
      this->vbat = 0.0f;
      this->actual_imax = 0.0f;
      this->tor_rc = 0l;
      this->water_ingress = 0l;
      this->pwr_relay = 0l;
    }
  }

  // field types and members
  using _header_type =
    std_msgs::msg::Header_<ContainerAllocator>;
  _header_type header;
  using _status_type =
    int32_t;
  _status_type status;
  using _temp_type =
    float;
  _temp_type temp;
  using _actual_pwm_type =
    int32_t;
  _actual_pwm_type actual_pwm;
  using _actual_pwm_rc_type =
    int32_t;
  _actual_pwm_rc_type actual_pwm_rc;
  using _pwm_source_type =
    int32_t;
  _pwm_source_type pwm_source;
  using _ibat_type =
    float;
  _ibat_type ibat;
  using _imot_type =
    float;
  _imot_type imot;
  using _vbat_type =
    float;
  _vbat_type vbat;
  using _actual_imax_type =
    float;
  _actual_imax_type actual_imax;
  using _tor_rc_type =
    int32_t;
  _tor_rc_type tor_rc;
  using _water_ingress_type =
    int32_t;
  _water_ingress_type water_ingress;
  using _pwr_relay_type =
    int32_t;
  _pwr_relay_type pwr_relay;

  // setters for named parameter idiom
  Type & set__header(
    const std_msgs::msg::Header_<ContainerAllocator> & _arg)
  {
    this->header = _arg;
    return *this;
  }
  Type & set__status(
    const int32_t & _arg)
  {
    this->status = _arg;
    return *this;
  }
  Type & set__temp(
    const float & _arg)
  {
    this->temp = _arg;
    return *this;
  }
  Type & set__actual_pwm(
    const int32_t & _arg)
  {
    this->actual_pwm = _arg;
    return *this;
  }
  Type & set__actual_pwm_rc(
    const int32_t & _arg)
  {
    this->actual_pwm_rc = _arg;
    return *this;
  }
  Type & set__pwm_source(
    const int32_t & _arg)
  {
    this->pwm_source = _arg;
    return *this;
  }
  Type & set__ibat(
    const float & _arg)
  {
    this->ibat = _arg;
    return *this;
  }
  Type & set__imot(
    const float & _arg)
  {
    this->imot = _arg;
    return *this;
  }
  Type & set__vbat(
    const float & _arg)
  {
    this->vbat = _arg;
    return *this;
  }
  Type & set__actual_imax(
    const float & _arg)
  {
    this->actual_imax = _arg;
    return *this;
  }
  Type & set__tor_rc(
    const int32_t & _arg)
  {
    this->tor_rc = _arg;
    return *this;
  }
  Type & set__water_ingress(
    const int32_t & _arg)
  {
    this->water_ingress = _arg;
    return *this;
  }
  Type & set__pwr_relay(
    const int32_t & _arg)
  {
    this->pwr_relay = _arg;
    return *this;
  }

  // constant declarations

  // pointer types
  using RawPtr =
    catarob_interfaces::msg::HullStatus_<ContainerAllocator> *;
  using ConstRawPtr =
    const catarob_interfaces::msg::HullStatus_<ContainerAllocator> *;
  using SharedPtr =
    std::shared_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator>>;
  using ConstSharedPtr =
    std::shared_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator> const>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::msg::HullStatus_<ContainerAllocator>>>
  using UniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator>, Deleter>;

  using UniquePtr = UniquePtrWithDeleter<>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::msg::HullStatus_<ContainerAllocator>>>
  using ConstUniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator> const, Deleter>;
  using ConstUniquePtr = ConstUniquePtrWithDeleter<>;

  using WeakPtr =
    std::weak_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator>>;
  using ConstWeakPtr =
    std::weak_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator> const>;

  // pointer types similar to ROS 1, use SharedPtr / ConstSharedPtr instead
  // NOTE: Can't use 'using' here because GNU C++ can't parse attributes properly
  typedef DEPRECATED__catarob_interfaces__msg__HullStatus
    std::shared_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator>>
    Ptr;
  typedef DEPRECATED__catarob_interfaces__msg__HullStatus
    std::shared_ptr<catarob_interfaces::msg::HullStatus_<ContainerAllocator> const>
    ConstPtr;

  // comparison operators
  bool operator==(const HullStatus_ & other) const
  {
    if (this->header != other.header) {
      return false;
    }
    if (this->status != other.status) {
      return false;
    }
    if (this->temp != other.temp) {
      return false;
    }
    if (this->actual_pwm != other.actual_pwm) {
      return false;
    }
    if (this->actual_pwm_rc != other.actual_pwm_rc) {
      return false;
    }
    if (this->pwm_source != other.pwm_source) {
      return false;
    }
    if (this->ibat != other.ibat) {
      return false;
    }
    if (this->imot != other.imot) {
      return false;
    }
    if (this->vbat != other.vbat) {
      return false;
    }
    if (this->actual_imax != other.actual_imax) {
      return false;
    }
    if (this->tor_rc != other.tor_rc) {
      return false;
    }
    if (this->water_ingress != other.water_ingress) {
      return false;
    }
    if (this->pwr_relay != other.pwr_relay) {
      return false;
    }
    return true;
  }
  bool operator!=(const HullStatus_ & other) const
  {
    return !this->operator==(other);
  }
};  // struct HullStatus_

// alias to use template instance with default allocator
using HullStatus =
  catarob_interfaces::msg::HullStatus_<std::allocator<void>>;

// constant definitions

}  // namespace msg

}  // namespace catarob_interfaces

#endif  // CATAROB_INTERFACES__MSG__DETAIL__HULL_STATUS__STRUCT_HPP_

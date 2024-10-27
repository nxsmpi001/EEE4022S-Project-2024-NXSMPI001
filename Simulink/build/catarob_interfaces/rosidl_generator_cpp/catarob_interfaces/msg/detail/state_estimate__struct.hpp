// generated from rosidl_generator_cpp/resource/idl__struct.hpp.em
// with input from catarob_interfaces:msg\StateEstimate.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__STRUCT_HPP_
#define CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__STRUCT_HPP_

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "rosidl_runtime_cpp/bounded_vector.hpp"
#include "rosidl_runtime_cpp/message_initialization.hpp"


#ifndef _WIN32
# define DEPRECATED__catarob_interfaces__msg__StateEstimate __attribute__((deprecated))
#else
# define DEPRECATED__catarob_interfaces__msg__StateEstimate __declspec(deprecated)
#endif

namespace catarob_interfaces
{

namespace msg
{

// message struct
template<class ContainerAllocator>
struct StateEstimate_
{
  using Type = StateEstimate_<ContainerAllocator>;

  explicit StateEstimate_(rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  {
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->x = 0.0;
      this->y = 0.0;
      this->psi = 0.0;
      this->u = 0.0;
      this->v = 0.0;
      this->r = 0.0;
      std::fill<typename std::array<double, 36>::iterator, double>(this->covariance.begin(), this->covariance.end(), 0.0);
    }
  }

  explicit StateEstimate_(const ContainerAllocator & _alloc, rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  : covariance(_alloc)
  {
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->x = 0.0;
      this->y = 0.0;
      this->psi = 0.0;
      this->u = 0.0;
      this->v = 0.0;
      this->r = 0.0;
      std::fill<typename std::array<double, 36>::iterator, double>(this->covariance.begin(), this->covariance.end(), 0.0);
    }
  }

  // field types and members
  using _x_type =
    double;
  _x_type x;
  using _y_type =
    double;
  _y_type y;
  using _psi_type =
    double;
  _psi_type psi;
  using _u_type =
    double;
  _u_type u;
  using _v_type =
    double;
  _v_type v;
  using _r_type =
    double;
  _r_type r;
  using _covariance_type =
    std::array<double, 36>;
  _covariance_type covariance;

  // setters for named parameter idiom
  Type & set__x(
    const double & _arg)
  {
    this->x = _arg;
    return *this;
  }
  Type & set__y(
    const double & _arg)
  {
    this->y = _arg;
    return *this;
  }
  Type & set__psi(
    const double & _arg)
  {
    this->psi = _arg;
    return *this;
  }
  Type & set__u(
    const double & _arg)
  {
    this->u = _arg;
    return *this;
  }
  Type & set__v(
    const double & _arg)
  {
    this->v = _arg;
    return *this;
  }
  Type & set__r(
    const double & _arg)
  {
    this->r = _arg;
    return *this;
  }
  Type & set__covariance(
    const std::array<double, 36> & _arg)
  {
    this->covariance = _arg;
    return *this;
  }

  // constant declarations

  // pointer types
  using RawPtr =
    catarob_interfaces::msg::StateEstimate_<ContainerAllocator> *;
  using ConstRawPtr =
    const catarob_interfaces::msg::StateEstimate_<ContainerAllocator> *;
  using SharedPtr =
    std::shared_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator>>;
  using ConstSharedPtr =
    std::shared_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator> const>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::msg::StateEstimate_<ContainerAllocator>>>
  using UniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator>, Deleter>;

  using UniquePtr = UniquePtrWithDeleter<>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::msg::StateEstimate_<ContainerAllocator>>>
  using ConstUniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator> const, Deleter>;
  using ConstUniquePtr = ConstUniquePtrWithDeleter<>;

  using WeakPtr =
    std::weak_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator>>;
  using ConstWeakPtr =
    std::weak_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator> const>;

  // pointer types similar to ROS 1, use SharedPtr / ConstSharedPtr instead
  // NOTE: Can't use 'using' here because GNU C++ can't parse attributes properly
  typedef DEPRECATED__catarob_interfaces__msg__StateEstimate
    std::shared_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator>>
    Ptr;
  typedef DEPRECATED__catarob_interfaces__msg__StateEstimate
    std::shared_ptr<catarob_interfaces::msg::StateEstimate_<ContainerAllocator> const>
    ConstPtr;

  // comparison operators
  bool operator==(const StateEstimate_ & other) const
  {
    if (this->x != other.x) {
      return false;
    }
    if (this->y != other.y) {
      return false;
    }
    if (this->psi != other.psi) {
      return false;
    }
    if (this->u != other.u) {
      return false;
    }
    if (this->v != other.v) {
      return false;
    }
    if (this->r != other.r) {
      return false;
    }
    if (this->covariance != other.covariance) {
      return false;
    }
    return true;
  }
  bool operator!=(const StateEstimate_ & other) const
  {
    return !this->operator==(other);
  }
};  // struct StateEstimate_

// alias to use template instance with default allocator
using StateEstimate =
  catarob_interfaces::msg::StateEstimate_<std::allocator<void>>;

// constant definitions

}  // namespace msg

}  // namespace catarob_interfaces

#endif  // CATAROB_INTERFACES__MSG__DETAIL__STATE_ESTIMATE__STRUCT_HPP_

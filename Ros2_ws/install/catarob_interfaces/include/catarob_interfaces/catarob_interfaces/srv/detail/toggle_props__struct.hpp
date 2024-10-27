// generated from rosidl_generator_cpp/resource/idl__struct.hpp.em
// with input from catarob_interfaces:srv/ToggleProps.idl
// generated code does not contain a copyright notice

#ifndef CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__STRUCT_HPP_
#define CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__STRUCT_HPP_

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "rosidl_runtime_cpp/bounded_vector.hpp"
#include "rosidl_runtime_cpp/message_initialization.hpp"


#ifndef _WIN32
# define DEPRECATED__catarob_interfaces__srv__ToggleProps_Request __attribute__((deprecated))
#else
# define DEPRECATED__catarob_interfaces__srv__ToggleProps_Request __declspec(deprecated)
#endif

namespace catarob_interfaces
{

namespace srv
{

// message struct
template<class ContainerAllocator>
struct ToggleProps_Request_
{
  using Type = ToggleProps_Request_<ContainerAllocator>;

  explicit ToggleProps_Request_(rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  {
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->flag = false;
    }
  }

  explicit ToggleProps_Request_(const ContainerAllocator & _alloc, rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  {
    (void)_alloc;
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->flag = false;
    }
  }

  // field types and members
  using _flag_type =
    bool;
  _flag_type flag;

  // setters for named parameter idiom
  Type & set__flag(
    const bool & _arg)
  {
    this->flag = _arg;
    return *this;
  }

  // constant declarations

  // pointer types
  using RawPtr =
    catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator> *;
  using ConstRawPtr =
    const catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator> *;
  using SharedPtr =
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator>>;
  using ConstSharedPtr =
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator> const>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator>>>
  using UniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator>, Deleter>;

  using UniquePtr = UniquePtrWithDeleter<>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator>>>
  using ConstUniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator> const, Deleter>;
  using ConstUniquePtr = ConstUniquePtrWithDeleter<>;

  using WeakPtr =
    std::weak_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator>>;
  using ConstWeakPtr =
    std::weak_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator> const>;

  // pointer types similar to ROS 1, use SharedPtr / ConstSharedPtr instead
  // NOTE: Can't use 'using' here because GNU C++ can't parse attributes properly
  typedef DEPRECATED__catarob_interfaces__srv__ToggleProps_Request
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator>>
    Ptr;
  typedef DEPRECATED__catarob_interfaces__srv__ToggleProps_Request
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Request_<ContainerAllocator> const>
    ConstPtr;

  // comparison operators
  bool operator==(const ToggleProps_Request_ & other) const
  {
    if (this->flag != other.flag) {
      return false;
    }
    return true;
  }
  bool operator!=(const ToggleProps_Request_ & other) const
  {
    return !this->operator==(other);
  }
};  // struct ToggleProps_Request_

// alias to use template instance with default allocator
using ToggleProps_Request =
  catarob_interfaces::srv::ToggleProps_Request_<std::allocator<void>>;

// constant definitions

}  // namespace srv

}  // namespace catarob_interfaces


#ifndef _WIN32
# define DEPRECATED__catarob_interfaces__srv__ToggleProps_Response __attribute__((deprecated))
#else
# define DEPRECATED__catarob_interfaces__srv__ToggleProps_Response __declspec(deprecated)
#endif

namespace catarob_interfaces
{

namespace srv
{

// message struct
template<class ContainerAllocator>
struct ToggleProps_Response_
{
  using Type = ToggleProps_Response_<ContainerAllocator>;

  explicit ToggleProps_Response_(rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  {
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->state = false;
    }
  }

  explicit ToggleProps_Response_(const ContainerAllocator & _alloc, rosidl_runtime_cpp::MessageInitialization _init = rosidl_runtime_cpp::MessageInitialization::ALL)
  {
    (void)_alloc;
    if (rosidl_runtime_cpp::MessageInitialization::ALL == _init ||
      rosidl_runtime_cpp::MessageInitialization::ZERO == _init)
    {
      this->state = false;
    }
  }

  // field types and members
  using _state_type =
    bool;
  _state_type state;

  // setters for named parameter idiom
  Type & set__state(
    const bool & _arg)
  {
    this->state = _arg;
    return *this;
  }

  // constant declarations

  // pointer types
  using RawPtr =
    catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator> *;
  using ConstRawPtr =
    const catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator> *;
  using SharedPtr =
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator>>;
  using ConstSharedPtr =
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator> const>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator>>>
  using UniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator>, Deleter>;

  using UniquePtr = UniquePtrWithDeleter<>;

  template<typename Deleter = std::default_delete<
      catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator>>>
  using ConstUniquePtrWithDeleter =
    std::unique_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator> const, Deleter>;
  using ConstUniquePtr = ConstUniquePtrWithDeleter<>;

  using WeakPtr =
    std::weak_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator>>;
  using ConstWeakPtr =
    std::weak_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator> const>;

  // pointer types similar to ROS 1, use SharedPtr / ConstSharedPtr instead
  // NOTE: Can't use 'using' here because GNU C++ can't parse attributes properly
  typedef DEPRECATED__catarob_interfaces__srv__ToggleProps_Response
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator>>
    Ptr;
  typedef DEPRECATED__catarob_interfaces__srv__ToggleProps_Response
    std::shared_ptr<catarob_interfaces::srv::ToggleProps_Response_<ContainerAllocator> const>
    ConstPtr;

  // comparison operators
  bool operator==(const ToggleProps_Response_ & other) const
  {
    if (this->state != other.state) {
      return false;
    }
    return true;
  }
  bool operator!=(const ToggleProps_Response_ & other) const
  {
    return !this->operator==(other);
  }
};  // struct ToggleProps_Response_

// alias to use template instance with default allocator
using ToggleProps_Response =
  catarob_interfaces::srv::ToggleProps_Response_<std::allocator<void>>;

// constant definitions

}  // namespace srv

}  // namespace catarob_interfaces

namespace catarob_interfaces
{

namespace srv
{

struct ToggleProps
{
  using Request = catarob_interfaces::srv::ToggleProps_Request;
  using Response = catarob_interfaces::srv::ToggleProps_Response;
};

}  // namespace srv

}  // namespace catarob_interfaces

#endif  // CATAROB_INTERFACES__SRV__DETAIL__TOGGLE_PROPS__STRUCT_HPP_

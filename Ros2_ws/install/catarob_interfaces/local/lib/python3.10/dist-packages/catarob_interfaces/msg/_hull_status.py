# generated from rosidl_generator_py/resource/_idl.py.em
# with input from catarob_interfaces:msg/HullStatus.idl
# generated code does not contain a copyright notice


# Import statements for member types

import builtins  # noqa: E402, I100

import math  # noqa: E402, I100

import rosidl_parser.definition  # noqa: E402, I100


class Metaclass_HullStatus(type):
    """Metaclass of message 'HullStatus'."""

    _CREATE_ROS_MESSAGE = None
    _CONVERT_FROM_PY = None
    _CONVERT_TO_PY = None
    _DESTROY_ROS_MESSAGE = None
    _TYPE_SUPPORT = None

    __constants = {
    }

    @classmethod
    def __import_type_support__(cls):
        try:
            from rosidl_generator_py import import_type_support
            module = import_type_support('catarob_interfaces')
        except ImportError:
            import logging
            import traceback
            logger = logging.getLogger(
                'catarob_interfaces.msg.HullStatus')
            logger.debug(
                'Failed to import needed modules for type support:\n' +
                traceback.format_exc())
        else:
            cls._CREATE_ROS_MESSAGE = module.create_ros_message_msg__msg__hull_status
            cls._CONVERT_FROM_PY = module.convert_from_py_msg__msg__hull_status
            cls._CONVERT_TO_PY = module.convert_to_py_msg__msg__hull_status
            cls._TYPE_SUPPORT = module.type_support_msg__msg__hull_status
            cls._DESTROY_ROS_MESSAGE = module.destroy_ros_message_msg__msg__hull_status

            from std_msgs.msg import Header
            if Header.__class__._TYPE_SUPPORT is None:
                Header.__class__.__import_type_support__()

    @classmethod
    def __prepare__(cls, name, bases, **kwargs):
        # list constant names here so that they appear in the help text of
        # the message class under "Data and other attributes defined here:"
        # as well as populate each message instance
        return {
        }


class HullStatus(metaclass=Metaclass_HullStatus):
    """Message class 'HullStatus'."""

    __slots__ = [
        '_header',
        '_status',
        '_temp',
        '_actual_pwm',
        '_actual_pwm_rc',
        '_pwm_source',
        '_ibat',
        '_imot',
        '_vbat',
        '_actual_imax',
        '_tor_rc',
        '_water_ingress',
        '_pwr_relay',
    ]

    _fields_and_field_types = {
        'header': 'std_msgs/Header',
        'status': 'int32',
        'temp': 'float',
        'actual_pwm': 'int32',
        'actual_pwm_rc': 'int32',
        'pwm_source': 'int32',
        'ibat': 'float',
        'imot': 'float',
        'vbat': 'float',
        'actual_imax': 'float',
        'tor_rc': 'int32',
        'water_ingress': 'int32',
        'pwr_relay': 'int32',
    }

    SLOT_TYPES = (
        rosidl_parser.definition.NamespacedType(['std_msgs', 'msg'], 'Header'),  # noqa: E501
        rosidl_parser.definition.BasicType('int32'),  # noqa: E501
        rosidl_parser.definition.BasicType('float'),  # noqa: E501
        rosidl_parser.definition.BasicType('int32'),  # noqa: E501
        rosidl_parser.definition.BasicType('int32'),  # noqa: E501
        rosidl_parser.definition.BasicType('int32'),  # noqa: E501
        rosidl_parser.definition.BasicType('float'),  # noqa: E501
        rosidl_parser.definition.BasicType('float'),  # noqa: E501
        rosidl_parser.definition.BasicType('float'),  # noqa: E501
        rosidl_parser.definition.BasicType('float'),  # noqa: E501
        rosidl_parser.definition.BasicType('int32'),  # noqa: E501
        rosidl_parser.definition.BasicType('int32'),  # noqa: E501
        rosidl_parser.definition.BasicType('int32'),  # noqa: E501
    )

    def __init__(self, **kwargs):
        assert all('_' + key in self.__slots__ for key in kwargs.keys()), \
            'Invalid arguments passed to constructor: %s' % \
            ', '.join(sorted(k for k in kwargs.keys() if '_' + k not in self.__slots__))
        from std_msgs.msg import Header
        self.header = kwargs.get('header', Header())
        self.status = kwargs.get('status', int())
        self.temp = kwargs.get('temp', float())
        self.actual_pwm = kwargs.get('actual_pwm', int())
        self.actual_pwm_rc = kwargs.get('actual_pwm_rc', int())
        self.pwm_source = kwargs.get('pwm_source', int())
        self.ibat = kwargs.get('ibat', float())
        self.imot = kwargs.get('imot', float())
        self.vbat = kwargs.get('vbat', float())
        self.actual_imax = kwargs.get('actual_imax', float())
        self.tor_rc = kwargs.get('tor_rc', int())
        self.water_ingress = kwargs.get('water_ingress', int())
        self.pwr_relay = kwargs.get('pwr_relay', int())

    def __repr__(self):
        typename = self.__class__.__module__.split('.')
        typename.pop()
        typename.append(self.__class__.__name__)
        args = []
        for s, t in zip(self.__slots__, self.SLOT_TYPES):
            field = getattr(self, s)
            fieldstr = repr(field)
            # We use Python array type for fields that can be directly stored
            # in them, and "normal" sequences for everything else.  If it is
            # a type that we store in an array, strip off the 'array' portion.
            if (
                isinstance(t, rosidl_parser.definition.AbstractSequence) and
                isinstance(t.value_type, rosidl_parser.definition.BasicType) and
                t.value_type.typename in ['float', 'double', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64']
            ):
                if len(field) == 0:
                    fieldstr = '[]'
                else:
                    assert fieldstr.startswith('array(')
                    prefix = "array('X', "
                    suffix = ')'
                    fieldstr = fieldstr[len(prefix):-len(suffix)]
            args.append(s[1:] + '=' + fieldstr)
        return '%s(%s)' % ('.'.join(typename), ', '.join(args))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self.header != other.header:
            return False
        if self.status != other.status:
            return False
        if self.temp != other.temp:
            return False
        if self.actual_pwm != other.actual_pwm:
            return False
        if self.actual_pwm_rc != other.actual_pwm_rc:
            return False
        if self.pwm_source != other.pwm_source:
            return False
        if self.ibat != other.ibat:
            return False
        if self.imot != other.imot:
            return False
        if self.vbat != other.vbat:
            return False
        if self.actual_imax != other.actual_imax:
            return False
        if self.tor_rc != other.tor_rc:
            return False
        if self.water_ingress != other.water_ingress:
            return False
        if self.pwr_relay != other.pwr_relay:
            return False
        return True

    @classmethod
    def get_fields_and_field_types(cls):
        from copy import copy
        return copy(cls._fields_and_field_types)

    @builtins.property
    def header(self):
        """Message field 'header'."""
        return self._header

    @header.setter
    def header(self, value):
        if __debug__:
            from std_msgs.msg import Header
            assert \
                isinstance(value, Header), \
                "The 'header' field must be a sub message of type 'Header'"
        self._header = value

    @builtins.property
    def status(self):
        """Message field 'status'."""
        return self._status

    @status.setter
    def status(self, value):
        if __debug__:
            assert \
                isinstance(value, int), \
                "The 'status' field must be of type 'int'"
            assert value >= -2147483648 and value < 2147483648, \
                "The 'status' field must be an integer in [-2147483648, 2147483647]"
        self._status = value

    @builtins.property
    def temp(self):
        """Message field 'temp'."""
        return self._temp

    @temp.setter
    def temp(self, value):
        if __debug__:
            assert \
                isinstance(value, float), \
                "The 'temp' field must be of type 'float'"
            assert not (value < -3.402823466e+38 or value > 3.402823466e+38) or math.isinf(value), \
                "The 'temp' field must be a float in [-3.402823466e+38, 3.402823466e+38]"
        self._temp = value

    @builtins.property
    def actual_pwm(self):
        """Message field 'actual_pwm'."""
        return self._actual_pwm

    @actual_pwm.setter
    def actual_pwm(self, value):
        if __debug__:
            assert \
                isinstance(value, int), \
                "The 'actual_pwm' field must be of type 'int'"
            assert value >= -2147483648 and value < 2147483648, \
                "The 'actual_pwm' field must be an integer in [-2147483648, 2147483647]"
        self._actual_pwm = value

    @builtins.property
    def actual_pwm_rc(self):
        """Message field 'actual_pwm_rc'."""
        return self._actual_pwm_rc

    @actual_pwm_rc.setter
    def actual_pwm_rc(self, value):
        if __debug__:
            assert \
                isinstance(value, int), \
                "The 'actual_pwm_rc' field must be of type 'int'"
            assert value >= -2147483648 and value < 2147483648, \
                "The 'actual_pwm_rc' field must be an integer in [-2147483648, 2147483647]"
        self._actual_pwm_rc = value

    @builtins.property
    def pwm_source(self):
        """Message field 'pwm_source'."""
        return self._pwm_source

    @pwm_source.setter
    def pwm_source(self, value):
        if __debug__:
            assert \
                isinstance(value, int), \
                "The 'pwm_source' field must be of type 'int'"
            assert value >= -2147483648 and value < 2147483648, \
                "The 'pwm_source' field must be an integer in [-2147483648, 2147483647]"
        self._pwm_source = value

    @builtins.property
    def ibat(self):
        """Message field 'ibat'."""
        return self._ibat

    @ibat.setter
    def ibat(self, value):
        if __debug__:
            assert \
                isinstance(value, float), \
                "The 'ibat' field must be of type 'float'"
            assert not (value < -3.402823466e+38 or value > 3.402823466e+38) or math.isinf(value), \
                "The 'ibat' field must be a float in [-3.402823466e+38, 3.402823466e+38]"
        self._ibat = value

    @builtins.property
    def imot(self):
        """Message field 'imot'."""
        return self._imot

    @imot.setter
    def imot(self, value):
        if __debug__:
            assert \
                isinstance(value, float), \
                "The 'imot' field must be of type 'float'"
            assert not (value < -3.402823466e+38 or value > 3.402823466e+38) or math.isinf(value), \
                "The 'imot' field must be a float in [-3.402823466e+38, 3.402823466e+38]"
        self._imot = value

    @builtins.property
    def vbat(self):
        """Message field 'vbat'."""
        return self._vbat

    @vbat.setter
    def vbat(self, value):
        if __debug__:
            assert \
                isinstance(value, float), \
                "The 'vbat' field must be of type 'float'"
            assert not (value < -3.402823466e+38 or value > 3.402823466e+38) or math.isinf(value), \
                "The 'vbat' field must be a float in [-3.402823466e+38, 3.402823466e+38]"
        self._vbat = value

    @builtins.property
    def actual_imax(self):
        """Message field 'actual_imax'."""
        return self._actual_imax

    @actual_imax.setter
    def actual_imax(self, value):
        if __debug__:
            assert \
                isinstance(value, float), \
                "The 'actual_imax' field must be of type 'float'"
            assert not (value < -3.402823466e+38 or value > 3.402823466e+38) or math.isinf(value), \
                "The 'actual_imax' field must be a float in [-3.402823466e+38, 3.402823466e+38]"
        self._actual_imax = value

    @builtins.property
    def tor_rc(self):
        """Message field 'tor_rc'."""
        return self._tor_rc

    @tor_rc.setter
    def tor_rc(self, value):
        if __debug__:
            assert \
                isinstance(value, int), \
                "The 'tor_rc' field must be of type 'int'"
            assert value >= -2147483648 and value < 2147483648, \
                "The 'tor_rc' field must be an integer in [-2147483648, 2147483647]"
        self._tor_rc = value

    @builtins.property
    def water_ingress(self):
        """Message field 'water_ingress'."""
        return self._water_ingress

    @water_ingress.setter
    def water_ingress(self, value):
        if __debug__:
            assert \
                isinstance(value, int), \
                "The 'water_ingress' field must be of type 'int'"
            assert value >= -2147483648 and value < 2147483648, \
                "The 'water_ingress' field must be an integer in [-2147483648, 2147483647]"
        self._water_ingress = value

    @builtins.property
    def pwr_relay(self):
        """Message field 'pwr_relay'."""
        return self._pwr_relay

    @pwr_relay.setter
    def pwr_relay(self, value):
        if __debug__:
            assert \
                isinstance(value, int), \
                "The 'pwr_relay' field must be of type 'int'"
            assert value >= -2147483648 and value < 2147483648, \
                "The 'pwr_relay' field must be an integer in [-2147483648, 2147483647]"
        self._pwr_relay = value

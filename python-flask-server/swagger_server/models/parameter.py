# coding: utf-8

from __future__ import absolute_import
from datetime import date, datetime  # noqa: F401

from typing import List, Dict  # noqa: F401

from swagger_server.models.base_model_ import Model
from swagger_server import util


class Parameter(Model):
    """NOTE: This class is auto generated by the swagger code generator program.

    Do not edit the class manually.
    """

    def __init__(self, name: str=None, type: str=None, default: str=None, biolink_class: str=None, allowed_values: List[str]=None, allowed_range: List[float]=None, suggested_values: str=None, lookup_url: str=None):  # noqa: E501
        """Parameter - a model defined in Swagger

        :param name: The name of this Parameter.  # noqa: E501
        :type name: str
        :param type: The type of this Parameter.  # noqa: E501
        :type type: str
        :param default: The default of this Parameter.  # noqa: E501
        :type default: str
        :param biolink_class: The biolink_class of this Parameter.  # noqa: E501
        :type biolink_class: str
        :param allowed_values: The allowed_values of this Parameter.  # noqa: E501
        :type allowed_values: List[str]
        :param allowed_range: The allowed_range of this Parameter.  # noqa: E501
        :type allowed_range: List[float]
        :param suggested_values: The suggested_values of this Parameter.  # noqa: E501
        :type suggested_values: str
        :param lookup_url: The lookup_url of this Parameter.  # noqa: E501
        :type lookup_url: str
        """
        self.swagger_types = {
            'name': str,
            'type': str,
            'default': str,
            'biolink_class': str,
            'allowed_values': List[str],
            'allowed_range': List[float],
            'suggested_values': str,
            'lookup_url': str
        }

        self.attribute_map = {
            'name': 'name',
            'type': 'type',
            'default': 'default',
            'biolink_class': 'biolink_class',
            'allowed_values': 'allowed_values',
            'allowed_range': 'allowed_range',
            'suggested_values': 'suggested_values',
            'lookup_url': 'lookup_url'
        }

        self._name = name
        self._type = type
        self._default = default
        self._biolink_class = biolink_class
        self._allowed_values = allowed_values
        self._allowed_range = allowed_range
        self._suggested_values = suggested_values
        self._lookup_url = lookup_url

    @classmethod
    def from_dict(cls, dikt) -> 'Parameter':
        """Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The parameter of this Parameter.  # noqa: E501
        :rtype: Parameter
        """
        return util.deserialize_model(dikt, cls)

    @property
    def name(self) -> str:
        """Gets the name of this Parameter.

        Name of the parameter.  # noqa: E501

        :return: The name of this Parameter.
        :rtype: str
        """
        return self._name

    @name.setter
    def name(self, name: str):
        """Sets the name of this Parameter.

        Name of the parameter.  # noqa: E501

        :param name: The name of this Parameter.
        :type name: str
        """
        if name is None:
            raise ValueError("Invalid value for `name`, must not be `None`")  # noqa: E501

        self._name = name

    @property
    def type(self) -> str:
        """Gets the type of this Parameter.

        Type of the parameter, one of 'Boolean', 'int', 'double', 'string'.  # noqa: E501

        :return: The type of this Parameter.
        :rtype: str
        """
        return self._type

    @type.setter
    def type(self, type: str):
        """Sets the type of this Parameter.

        Type of the parameter, one of 'Boolean', 'int', 'double', 'string'.  # noqa: E501

        :param type: The type of this Parameter.
        :type type: str
        """
        allowed_values = ["Boolean", "int", "double", "string"]  # noqa: E501
        if type not in allowed_values:
            raise ValueError(
                "Invalid value for `type` ({0}), must be one of {1}"
                .format(type, allowed_values)
            )

        self._type = type

    @property
    def default(self) -> str:
        """Gets the default of this Parameter.

        Default value of the parameter.  # noqa: E501

        :return: The default of this Parameter.
        :rtype: str
        """
        return self._default

    @default.setter
    def default(self, default: str):
        """Sets the default of this Parameter.

        Default value of the parameter.  # noqa: E501

        :param default: The default of this Parameter.
        :type default: str
        """
        if default is None:
            raise ValueError("Invalid value for `default`, must not be `None`")  # noqa: E501

        self._default = default

    @property
    def biolink_class(self) -> str:
        """Gets the biolink_class of this Parameter.

        Biolink class of the paramater. Applicable to producers only and only one parameter can have a biolink class.  # noqa: E501

        :return: The biolink_class of this Parameter.
        :rtype: str
        """
        return self._biolink_class

    @biolink_class.setter
    def biolink_class(self, biolink_class: str):
        """Sets the biolink_class of this Parameter.

        Biolink class of the paramater. Applicable to producers only and only one parameter can have a biolink class.  # noqa: E501

        :param biolink_class: The biolink_class of this Parameter.
        :type biolink_class: str
        """

        self._biolink_class = biolink_class

    @property
    def allowed_values(self) -> List[str]:
        """Gets the allowed_values of this Parameter.

        Allowed values for the parameter.  # noqa: E501

        :return: The allowed_values of this Parameter.
        :rtype: List[str]
        """
        return self._allowed_values

    @allowed_values.setter
    def allowed_values(self, allowed_values: List[str]):
        """Sets the allowed_values of this Parameter.

        Allowed values for the parameter.  # noqa: E501

        :param allowed_values: The allowed_values of this Parameter.
        :type allowed_values: List[str]
        """

        self._allowed_values = allowed_values

    @property
    def allowed_range(self) -> List[float]:
        """Gets the allowed_range of this Parameter.

        Allowed range for values of the parameter.  # noqa: E501

        :return: The allowed_range of this Parameter.
        :rtype: List[float]
        """
        return self._allowed_range

    @allowed_range.setter
    def allowed_range(self, allowed_range: List[float]):
        """Sets the allowed_range of this Parameter.

        Allowed range for values of the parameter.  # noqa: E501

        :param allowed_range: The allowed_range of this Parameter.
        :type allowed_range: List[float]
        """

        self._allowed_range = allowed_range

    @property
    def suggested_values(self) -> str:
        """Gets the suggested_values of this Parameter.

        Suggested value range for the parameter.  # noqa: E501

        :return: The suggested_values of this Parameter.
        :rtype: str
        """
        return self._suggested_values

    @suggested_values.setter
    def suggested_values(self, suggested_values: str):
        """Sets the suggested_values of this Parameter.

        Suggested value range for the parameter.  # noqa: E501

        :param suggested_values: The suggested_values of this Parameter.
        :type suggested_values: str
        """

        self._suggested_values = suggested_values

    @property
    def lookup_url(self) -> str:
        """Gets the lookup_url of this Parameter.

        URL to search for suitable parameter values.  # noqa: E501

        :return: The lookup_url of this Parameter.
        :rtype: str
        """
        return self._lookup_url

    @lookup_url.setter
    def lookup_url(self, lookup_url: str):
        """Sets the lookup_url of this Parameter.

        URL to search for suitable parameter values.  # noqa: E501

        :param lookup_url: The lookup_url of this Parameter.
        :type lookup_url: str
        """

        self._lookup_url = lookup_url

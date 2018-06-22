# -*- coding: utf-8 -*-

""" Custom exceptions """

class BadInputException(Exception):
    """ Python error if no suitable input are provided.
    """
    def __init__(self, message):
        super().__init__(message)

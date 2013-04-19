#!/usr/bin/python
# -*- coding: utf-8 -*-

import re


class InvalidInputError(Exception):

    pass


def first_mol(sdf):
    """given one sdf file or pseudo file or string, use the 'read' method to
    read only the first MOL section if sdf is a file or pseudo file, and get
    the first MOL section from input if input is a string"""

    try:
        content = sdf.read(160000)
    except:
        content = sdf
    m_end = re.compile(r'^M\s+END', re.MULTILINE)
    match = m_end.search(content)
    if not match:
        raise InvalidInputError
    return content[:match.end()] + '\n$$$$'



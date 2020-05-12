"""
Define custom exceptions for PneumoCaT2.
Carmen Sheppard 2020
"""
class CtvdbError(Exception):
    """create custom error for missing or mismatching information in DB. Allows custom messages"""

    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"CtvdbError: {self.message}"
        else:
            return "CtvdbError: check CTV.db and folder integrity, missing or mismatching " \
                   "information may be present.\n"

class CtvdbFileError(Exception):
    """create custom error for missing or mismatching files in ctvdb folder. Allows custom messages"""

    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"CtvdbFileError: {self.message}"
        else:
            return "CtvdbFileError: check folder integrity for missing or misnamed " \
                   " reference files.\n"
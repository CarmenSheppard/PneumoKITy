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
                   "information may be present."


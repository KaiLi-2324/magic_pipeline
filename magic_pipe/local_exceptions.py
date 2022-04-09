class SamplesNotFoundError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class FormatLineNotFoundError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class TempDirNotCleanError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class EmptyResultError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class UnExpectedError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class InvalidSampleTypeError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class GenesParsingError(BaseException):
    def __init__(self, genes):
        self.genes = genes

    def __str__(self):
        return "Unexpected Error occurred when parsing the following genes {}".format(",".join(self.genes))


class GenesPermNumberNotEqualError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class SamplesNotEqualError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class VersionError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info
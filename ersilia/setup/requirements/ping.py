from ... import throw_ersilia_exception


class PingRequirement(object):
    def __init__(self):
        self.is_connected()

    @throw_ersilia_exception
    def is_connected(self):
        pass
        # TODO

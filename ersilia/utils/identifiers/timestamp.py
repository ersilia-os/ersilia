from datetime import datetime


class TimeStampIdentifier(object):
    def __init__(self):
        self.stamp = datetime.now()

    def encode(self):
        return self.stamp.strftime("%Y%m%d%H%M%S")

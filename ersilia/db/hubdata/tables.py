class DynamoDbTable(object):
    """
    Base class for DynamoDB tables.
    """

    def __init__(self):
        pass


class PredictionsTable(DynamoDbTable):
    """
    Table for storing predictions.
    """

    def __init__(self):
        pass


class ModelsTable(DynamoDbTable):
    """
    Table for storing models.
    """

    def __init__(self):
        pass

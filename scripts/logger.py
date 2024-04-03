import logging

# Custom formatter for the logging module
COLORS = {
    'WARNING': '\033[93m',
    'INFO': '\033[92m',
    'DEBUG': '\033[94m',
    'CRITICAL': '\033[91m',
    'ERROR': '\033[91m',
    'ENDC': '\033[0m',
}


class CustomFormatter(logging.Formatter):
    def format(self, record):
        levelname = record.levelname
        message = logging.Formatter.format(self, record)
        return f"{COLORS.get(levelname, '')}{message}{COLORS['ENDC']}"


def custom_logger(name=__name__):
    """
    Creates a custom logger with a custom name. Custom 
    logger with a specific name. It sets a custom DEBUG and formatters that
    include the timestamp, log level, and message. Every log level has 
    it own color this increases the readability.
    All previous handlers logger are removed.

    Args:
        name (str): The name of the logger as string, mostly it is __name__.

    Returns:
        logger(Logger): The custom logger.

    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    formatter = CustomFormatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.handlers = []
    logger.addHandler(handler)
    return logger

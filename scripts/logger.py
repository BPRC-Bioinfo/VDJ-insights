import logging

ENVIRONMENT_LEVEL = 25
SUCCESS_LEVEL = 26
NOTICE_LEVEL = 27
TRACE_LEVEL = 5

logging.addLevelName(ENVIRONMENT_LEVEL, "ENVIRONMENT")
logging.addLevelName(SUCCESS_LEVEL, "SUCCESS")
logging.addLevelName(NOTICE_LEVEL, "NOTICE")
logging.addLevelName(TRACE_LEVEL, "TRACE")

COLORS = {
    'WARNING': '\033[93m',      # Yellow
    'INFO': '\033[92m',         # Green
    'DEBUG': '\033[94m',        # Blue
    'CRITICAL': '\033[91m',     # Red
    'ERROR': '\033[91m',        # Red
    'ENVIRONMENT': '\033[95m',  # Purple
    'SUCCESS': '\033[92m',      # Green
    'NOTICE': '\033[96m',       # Cyan
    'TRACE': '\033[90m',        # Grey
    'ENDC': '\033[0m',          # Reset
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
    its own color, which increases readability.
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
    formatter = CustomFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.handlers = []
    logger.addHandler(handler)

    logging.Logger.environment = lambda self, message, * \
        args, **kws: self._log(ENVIRONMENT_LEVEL, message, args, **kws)
    logging.Logger.success = lambda self, message, * \
        args, **kws: self._log(SUCCESS_LEVEL, message, args, **kws)
    logging.Logger.notice = lambda self, message, * \
        args, **kws: self._log(NOTICE_LEVEL, message, args, **kws)
    logging.Logger.trace = lambda self, message, * \
        args, **kws: self._log(TRACE_LEVEL, message, args, **kws)

    return logger

import logging
from pathlib import Path

ENVIRONMENT_LEVEL = 25
SUCCESS_LEVEL = 26
NOTICE_LEVEL = 27
TRACE_LEVEL = 5

logging.addLevelName(ENVIRONMENT_LEVEL, "ENVIRONMENT")
logging.addLevelName(SUCCESS_LEVEL, "SUCCESS")
logging.addLevelName(NOTICE_LEVEL, "NOTICE")
logging.addLevelName(TRACE_LEVEL, "TRACE")

COLORS = {
    'WARNING': '\033[93m',    # Yellow
    'INFO': '\033[92m',       # Green
    'DEBUG': '\033[94m',      # Blue
    'CRITICAL': '\033[91m',   # Red
    'ERROR': '\033[91m',      # Red
    'ENVIRONMENT': '\033[95m',  # Purple
    'SUCCESS': '\033[92m',    # Green
    'NOTICE': '\033[96m',     # Cyan
    'TRACE': '\033[90m',      # Grey
    'ENDC': '\033[0m',        # Reset
}


class CustomFormatter(logging.Formatter):
    def format(self, record):
        levelname = record.levelname
        message = logging.Formatter.format(self, record)
        return f"{COLORS.get(levelname, '')}{message}{COLORS['ENDC']}"


def console_logger(name=__name__):
    """
    Creates a logger that prints log messages to the console with color formatting.
    Args:
        name (str): The name of the logger as a string.
    Returns:
        logger (Logger): Logger configured to print to the console.
    """
    logger = logging.getLogger(f'{name}_console')
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_formatter = CustomFormatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    stream_handler.setFormatter(stream_formatter)

    logger.handlers = []
    logger.addHandler(stream_handler)

    _add_custom_log_levels(logger)

    return logger


def file_logger(name=__name__, log_file='logging.log'):
    """
    Creates a logger that writes log messages to a file.
    Args:
        name (str): The name of the logger as a string.
        log_file (str): The filename for the log file output.
    Returns:
        logger (Logger): Logger configured to write to a file.
    """
    cwd = Path.cwd()
    log_file_path = cwd / log_file

    logger = logging.getLogger(f'{name}_file')
    logger.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(file_formatter)

    logger.handlers = []
    logger.addHandler(file_handler)

    _add_custom_log_levels(logger)

    return logger


def _add_custom_log_levels(logger):
    """
    Adds custom log levels to the logger.
    Args:
        logger (Logger): Logger to which custom log levels will be added.
    """
    logging.Logger.environment = lambda self, message, *args, **kws: self._log(ENVIRONMENT_LEVEL, message, args, **kws)
    logging.Logger.success = lambda self, message, *args, **kws: self._log(SUCCESS_LEVEL, message, args, **kws)
    logging.Logger.notice = lambda self, message, *args, **kws: self._log(NOTICE_LEVEL, message, args, **kws)
    logging.Logger.trace = lambda self, message, *args, **kws: self._log(TRACE_LEVEL, message, args, **kws)


if __name__ == "__main__":
    console_log = console_logger(__name__)
    file_log = file_logger(__name__)

    file_log.warning("This is a warning message printed to the console")
    file_log.success("This is a success message printed to the console")
    file_log.notice("This is a notice message printed to the console")

    file_log.error("This is an error message written to the file")
    file_log.warning("This is a trace message written to the file")
    file_log.environment("This is an environment message written to the file")

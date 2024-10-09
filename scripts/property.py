import sys

from logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)

def log_error():
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if hasattr(e, 'value'):
                    file_log.error(f"Error '{func.__name__}': {e} | Input: {e.value}")
                else:
                    file_log.error(f"Error '{func.__name__}': {e}")
                sys.exit(1)
        return wrapper
    return decorator



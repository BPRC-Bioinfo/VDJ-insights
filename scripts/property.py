import sys

from logger import custom_logger

logger = custom_logger(__name__)

def log_error():
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if hasattr(e, 'value'):
                    logger.error(f"Error '{func.__name__}': {e} | Input: {e.value}")
                else:
                    logger.error(f"Error '{func.__name__}': {e}")
                sys.exit(1)
        return wrapper
    return decorator



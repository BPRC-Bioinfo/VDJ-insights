import sys

from logger import console_logger, file_logger

console_log = console_logger(__name__)
file_log = file_logger(__name__)

def log_error():
    """
    Decorator to log errors that occur in the decorated function.

    Returns:
        function: A decorator function that wraps the original function.
    """
    def decorator(func):
        """
        Inner decorator function.

        Args:
            func (function): The function to be decorated.

        Returns:
            function: The wrapped function with error logging.
        """
        def wrapper(*args, **kwargs):
            """
               Inner decorator function.

               Args:
                   func (function): The function to be decorated.

               Returns:
                   function: The wrapped function with error logging.
               """

            def wrapper(*args, **kwargs):
            """
            Wrapper function that executes the original function and logs any exceptions.

            Args:
                *args: Variable length argument list for the original function.
                **kwargs: Arbitrary keyword arguments for the original function.

            Returns:
                Any: The return value of the original function, if no exception occurs.

            Raises:
                SystemExit: Exits the program with status 1 if an exception occurs.
            """
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if hasattr(e, 'value'):
                    console_log.error(f"Error '{func.__name__}': {e} | Input: {e.value}")
                else:
                    console_log.error(f"Error '{func.__name__}': {e}")
                sys.exit(1)
        return wrapper
    return decorator



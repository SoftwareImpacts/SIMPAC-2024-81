from datetime import datetime
import functools
import inspect
import os
import warnings

# Global variable to track the last call signature
last_call_signature = None


def get_log_file_name(suffix: str = "log.txt") -> str:
    """
    Generate a log file name with a timestamp.

    Args:
        suffix (str, optional): The suffix to be added to the filename. Defaults to "log.txt".

    Returns:
        str: A string representing the log file name with a timestamp.
    """
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")
    log_file_name = f"{timestamp}_{suffix}"
    return log_file_name


def create_log_file(logdir: str = "", constant: bool = False) -> str:
    """
    Create a log file and return its absolute path.

    Args:
        logdir (str, optional): The directory where the log file will be created. Defaults to the current directory.
        constant (bool, optional): If True, the log file will have a constant name 'log.txt' instead of timestamped. Defaults to False.

    Returns:
        str: The absolute path to the created log file.
    """
    if not constant:
        logfile_name = logdir + get_log_file_name(suffix="log.txt")
    else:
        logfile_name = logdir + "log.txt"
    absolute_path = os.path.abspath(logfile_name)
    with open(absolute_path, "w") as log_file:
        log_file.write("")  # Clear the contents of the file by writing an empty string
    return absolute_path


def print_to_log(logfile: str, message: str, warning: bool = False) -> None:
    """
    Print a message to a log file with timestamp and optionally raise a warning.

    Args:
        logfile (str): The path to the log file where the message will be written.
        message (str): The message to be written to the log file.
        warning (bool, optional): If True, raise a warning with the message. Defaults to False.

    Returns:
        None
    """
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H:%M:%S")

    log_message = f"{timestamp} | {message}"

    with open(logfile, "a") as log_file:
        log_file.write(log_message + "\n")

    if warning:
        warnings.warn(log_message)

    print(log_message, flush=True)


def create_function_calls_log_file(
    logfile: str = "calls_log.txt", title: str = "**Program Initialization**"
):
    """
    Initializes or resets the log file used for function call logging. Optionally,
    a custom log file name and title can be specified.

    :param logfile: Name of the log file. Defaults to "calls_log.txt".
    :param title: Title to write at the top of the log file. Defaults to "**Program Initialization**".
    """
    try:
        os.remove(logfile)
    except:
        pass
    with open("calls_log.txt", "w") as log_file:
        log_file.write(title)


def track_call_depth(func):
    """
    Decorator that logs calls to the decorated function with timestamp and call signature.
    Differentiates consecutive calls to different functions, appending details to `calls_log.txt`.
    Utilizes `inspect` for call signature and `datetime` for timestamps.

    :param func: Function to decorate.
    :return: Wrapper function with logging capability.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global last_call_signature

        # Construct the current call signature
        frame = inspect.currentframe().f_back
        filename = os.path.basename(frame.f_code.co_filename).replace(".py", "")
        class_name = ""
        if "self" in frame.f_locals:
            class_name = frame.f_locals["self"].__class__.__name__ + "."
        current_call_signature = f"{filename}.{class_name}{func.__name__}"

        # Get the current timestamp
        now = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

        # Check if the current call signature matches the last one
        if current_call_signature != last_call_signature:
            # Construct the log message
            log_message = f"{now} | {current_call_signature}\n"

            # Open the log file and append the log message
            with open("calls_log.txt", "a") as log_file:
                log_file.write(log_message)

            last_call_signature = current_call_signature

        # Call the original function
        return func(*args, **kwargs)

    return wrapper

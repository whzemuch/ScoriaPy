import logging
from typing import Optional


def setup_logger(log_file=None, level=logging.INFO):
    """
    Create and configure a logger for Scanpy processing pipelines.

    This logger:
    - prints messages to the console with timestamps
    - optionally writes logs to a file
    - avoids duplicate handlers (important in Jupyter notebooks)
    - supports adjustable log levels (INFO, DEBUG, WARNING, etc.)

    Parameters
    ----------
    log_file : str or None, optional (default: None)
        If provided, log messages are also written to this file.
        If None, only console output is produced.
    level : int, optional (default: logging.INFO)
        Logging verbosity level. Typically one of:
        - logging.DEBUG
        - logging.INFO
        - logging.WARNING
        - logging.ERROR
        - logging.CRITICAL

    Returns
    -------
    logger : logging.Logger
        A configured logger instance that can be passed to functions
        (e.g., preprocessing pipelines) to produce timestamped logs.

    Notes
    -----
    Inside Jupyter notebooks, calling this function repeatedly will not
    produce duplicate log messages because the logger only adds handlers
    the first time it is constructed.

    Examples
    --------
    Basic usage:

    >>> logger = setup_logger()
    >>> logger.info("Starting preprocessing...")

    Logging to both console and file:

    >>> logger = setup_logger("scanpy_run.log")
    >>> logger.info("Running HVG selection...")

    Using inside a Scanpy pipeline:

    >>> logger = setup_logger("fibro_preprocess.log")
    >>> adata = run_scanpy_basic_pipeline(adata, logger=logger)
    >>> logger.info("Pipeline finished.")

    Changing verbosity:

    >>> logger = setup_logger(level=logging.DEBUG)
    >>> logger.debug("Debug-level message for troubleshooting.")
    """
    logger = logging.getLogger("scanpy_pipeline")
    logger.setLevel(level)

    # Avoid adding handlers multiple times (important for Jupyter)
    if logger.handlers:
        return logger

    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Optional file handler
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger
import logging

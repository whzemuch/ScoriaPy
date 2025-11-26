import logging
from typing import Optional


def setup_logger(log_file: Optional[str] = None, level=logging.INFO):
    """
    Create and configure the shared ``\"scoriapy\"`` logger.

    The logger is configured only once; subsequent calls reuse the same
    instance without adding duplicate handlers. By default, log messages are
    written to ``stderr``; if ``log_file`` is provided, they are also written
    to a file.

    Parameters
    ----------
    log_file
        Optional path to a log file. If given, a ``FileHandler`` is added in
        addition to the standard stream handler.
    level
        Logging level to use (e.g. ``logging.INFO``).

    Returns
    -------
    logging.Logger
        Configured logger instance named ``\"scoriapy\"``.
    """

    logger = logging.getLogger("scoriapy")
    logger.setLevel(level)

    if not logger.handlers:
        formatter = logging.Formatter("%(asctime)s — %(levelname)s — %(message)s")

        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        if log_file:
            fh = logging.FileHandler(log_file)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

    return logger

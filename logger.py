import logging
import os
import time

from config import LOG_LEVEL


class LogHandler(logging.StreamHandler):
    """ Получение и обработка лога."""

    def __init__(self):
        logging.StreamHandler.__init__(self)
        fmt = '%(asctime)s | %(filename)s | %(levelname)s | %(message)s'
        fmt_date = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(fmt, fmt_date)
        self.setFormatter(formatter)
        self.setLevel('DEBUG')


logs_dir = os.path.join(os.path.dirname(__file__))
log_file_path = os.path.join(logs_dir, f"{__name__}{time.time()}.log")

logger = logging.getLogger(__name__)
logger.setLevel(level=LOG_LEVEL)

handler = logging.FileHandler(log_file_path, mode='a')
logger.addHandler(handler)
logger.addHandler(LogHandler())

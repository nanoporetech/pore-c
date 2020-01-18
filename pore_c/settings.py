import logging
from logging.config import dictConfig


BASE_CONFIG = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {"simple": {"format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s"}},
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "level": "DEBUG",
            "formatter": "simple",
            "stream": "ext://sys.stderr",
        }
    },
    "loggers": {"pore_c": {"level": "ERROR", "handlers": ["console"], "propagate": False}},
    "root": {"level": "INFO", "handlers": ["console"]},
}


def setup_logging():
    dictConfig(BASE_CONFIG)
    return logging.getLogger("pore_c")

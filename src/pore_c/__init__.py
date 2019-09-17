import logging

from ._version import get_versions

logging.getLogger(__name__).addHandler(logging.NullHandler())

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

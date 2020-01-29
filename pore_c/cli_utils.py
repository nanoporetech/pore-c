import json
import re
from logging import getLogger
from pathlib import Path

import click


logger = getLogger(__name__)


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.
    """

    def list_commands(self, ctx):
        """List command names as they are in commands dict.
        """
        return self.commands.keys()


def command_line_json(ctx, param, value):
    # TODO: add support for json from file
    if value is None:
        return {}
    try:
        res = json.loads(value)
    except Exception as exc:  # noqa: F841
        logger.exception("Not valid json")
        raise
    return res


def expand_output_prefix(catalog):
    def inner(ctx, param, value):

        file_paths = catalog.generate_paths(value)
        ctx.meta["file_paths"] = file_paths
        return value

    return inner


def filename_matches_regex(regex_str: str):
    regex = re.compile(regex_str)

    def inner(ctx, param, value):

        value = Path(value)
        m = regex.match(value.name)
        if not m:
            raise click.BadParameter(f"Filename {value.name} doesn't match {regex}")
        return value

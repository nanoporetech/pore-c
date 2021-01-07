import json
import re
from logging import getLogger
from pathlib import Path

import click
from pysam import AlignmentFile


logger = getLogger(__name__)


class ExportDependentOption(click.Option):
    def __init__(self, *args, **kwds):
        self.export_formats: list = kwds.pop("export_formats")

        kwds["help"] = kwds.get("help", "") + "(required if export format is in {})".format(
            ",".join(self.export_formats)
        )
        super(ExportDependentOption, self).__init__(*args, **kwds)

    def handle_parse_result(self, ctx, opts, args):
        fmt = ctx.params.get("format", None)
        if fmt in self.export_formats:
            if self.name not in opts:
                raise click.UsageError("When exporting to {}, {} is a required parameter".format(fmt, self.name))
            else:
                self.prompt = None
        return super(ExportDependentOption, self).handle_parse_result(ctx, opts, args)


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


def pipeable_sam_input(ctx, param, value):
    if value == "-":
        if ctx.params["input_is_bam"]:
            mode = "rb"
        else:
            mode = "r"
        res = AlignmentFile(value, mode=mode)
    else:
        res = AlignmentFile(value)
    return res


def pipeable_sam_output(ctx, param, value):
    input_sam = ctx.params["input_sam"]

    if value.endswith(".bam") or (value == "-" and ctx.params["output_is_bam"]):
        mode = "wb"
    else:
        mode = "w"
    res = AlignmentFile(value, mode=mode, template=input_sam)
    return res

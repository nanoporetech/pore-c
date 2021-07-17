import sys

import click
import pandas as pd

from pore_c.cli import cli

ctx = click.Context(cli)


data = []


def _get_cmd_df(cmd, ctx):
    data = []
    for group_name, group in cli.commands.items():
        data.append({"group": group_name})
        for cmd_name, cmd in group.commands.items():
            data.append({"group": group_name, "command": cmd_name, "description": cmd.short_help, "_cmd": cmd})
    return pd.DataFrame(data)


def _create_toc(cmd_df, ctx):

    df = cmd_df.copy().fillna("")
    res = ["# Command Line Interface\n\n"]
    res.append("""```{contents} Table of Contents\n---\ndepth: 3\n---\n```\n\n""")
    for _, row in df.iterrows():
        _row = row.to_dict()
        _cmd = _row["_cmd"]
        if _cmd != "":
            res.append(f"### {row['command']}\n\n{row['description']}\n\n")
            usage = _cmd.get_usage(ctx).strip().replace("Usage:", "").strip()
            usage = f"pore_c {row['group']} {row['command']} {usage}"
            if _cmd.help:
                _help = _cmd.help
            else:
                _help = ""
            res.append(f"#### Usage\n```bash\n{usage}\n\n{_help}\n```\n")

            res.append("#### Parameters:\n\n")
            for param in _cmd.get_params(ctx):
                if param.param_type_name == "argument":
                    res.append(f"- *{param.name}* [required]\n")
                else:
                    opt, _help = param.get_help_record(ctx)
                    if opt != "--help":
                        res.append(f"- *{opt}*: {_help}\n")
            res.append("\n***\n")
        else:
            res.append(f"## {row['group']}\n")
    if res[-1] == "\n***\n":
        res = res[:-1]
    return "".join(res)


cmd_df = _get_cmd_df(cli, ctx)
toc = _create_toc(cmd_df, ctx)
sys.stdout.write(
    f"""

{toc}
"""
)

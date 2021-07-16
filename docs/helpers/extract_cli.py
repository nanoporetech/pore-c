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


def _create_toc(cmd_df):

    df = cmd_df[["group", "command", "description"]].copy()
    df.loc[~df["command"].isnull(), "group"] = ""

    return df.fillna("").set_index(["group"]).to_markdown()


def _get_cmd_table(cmd, ctx):
    for param in cmd.get_params(ctx):
        print(param.get_help_record(ctx))


# print(_get_cmd_table(cli, ctx))
# print(_get_cmd_df(cli, ctx).fillna("").drop("_cmd", 1).set_index(["group"]).to_markdown(tablefmt="github"))


cmd_df = _get_cmd_df(cli, ctx)

toc = _create_toc(cmd_df)


#
# for name, group in cli.commands.items():
#    print(name, group)
#    for cmd_name, cmd in group.commands.items():
#        data.append({"group": name, "command": cmd_name, "description": cmd.short_help})
#        print(cmd.short_help)
#        for param in cmd.get_params(ctx):
#            print(param.name, param.param_type_name)
#            print(param.get_help_record(ctx))
#            print(dir(param))
#

# df = pd.DataFrame(data).set_index(["group"])


with open(sys.argv[1], "w") as fh:

    fh.write(
        f"""

{toc}
    """
    )

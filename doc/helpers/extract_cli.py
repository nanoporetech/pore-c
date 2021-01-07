import sys

import pandas as pd

from pore_c.cli import cli


data = []
for name, group in cli.commands.items():
    for cmd_name, cmd in group.commands.items():
        data.append({"group": name, "command": cmd_name, "description": cmd.short_help})


df = pd.DataFrame(data).set_index(["group"])

with open(sys.argv[1], "w") as fh:
    df.to_markdown(fh)

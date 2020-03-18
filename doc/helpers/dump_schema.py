import pore_c.model as m
from pydantic.schema import schema
import json


top_level_schema = schema(
    [
    m.FragmentRecord,
    m.AlignmentRecord,
    m.PoreCRecord,
    m.PoreCContactRecord,
    m.PoreCConcatemerRecord
], title="PoreC Data Model")

print(json.dumps(top_level_schema, indent=2))

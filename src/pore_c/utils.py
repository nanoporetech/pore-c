import re

def kmg_bases_to_int(value: str) -> int:
    try:
        result = int(value)
    except:
        result = None
    if result is not None:
        return result


    value_re = re.compile(r"(\d+)([KkMmGg])b*")
    m = value_re.match(value)
    if not m:
        raise ValueError(f"Invalid string: {value}")
    value = int(m.group(1))
    exponent = {"k": 1e3, 'm': 1e6, 'g': 1e9}[m.group(2).lower()]
    return value * int(exponent)


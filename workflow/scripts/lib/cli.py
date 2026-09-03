"""Small helpers for constructing command arguments without shell fragments."""

from __future__ import annotations


def command_option_values(option: str, comma_separated_values: str) -> list[str]:
    """Expand comma-separated values into repeated, whitespace-trimmed options."""

    if not comma_separated_values:
        return []
    values = [
        value.strip()
        for value in comma_separated_values.split(",")
        if value.strip()
    ]
    return [token for value in values for token in (option, value)]

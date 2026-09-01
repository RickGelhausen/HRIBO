#!/usr/bin/env python
"""
Parsing and formatting of the GFF/GTF attribute column.

Kept dependency free -- standard library only -- because the scripts that need
it run in the mergetools environment, which has no biopython.

The scripts here previously each split the attribute column on both ";" and "="
at once and then paired the resulting items up two at a time. That works only
when every attribute has a non-empty value: an empty one such as "ORF_type=;"
collapses to a single item and shifts every following key onto the wrong value,
or runs off the end of the list entirely. Splitting on ";" first and only then on
the first "=" avoids that, and also survives a value that itself contains "=".

Author: Rick Gelhausen
"""

import re
import sys

# GTF2 writes attributes as: key "value"; key "value";
GTF2_PAIR = re.compile(r'^\s*(?P<key>\S+)\s+"(?P<value>[^"]*)"\s*$')

# GFF3 reserves capitalised attribute names for this defined set.  Prediction
# tools historically emitted names such as ``Prob`` and ``Evidence``; those are
# user attributes and must start with a lowercase character to remain valid
# GFF3.  Keep the spelling of the real reserved names while normalising every
# custom key.
GFF3_RESERVED_ATTRIBUTES = {
    key.lower(): key
    for key in (
        "ID",
        "Name",
        "Alias",
        "Parent",
        "Target",
        "Gap",
        "Derives_from",
        "Note",
        "Dbxref",
        "Ontology_term",
        "Is_circular",
    )
}


def split_attributes(attributes):
    """Ordered (key, value) pairs, preserving key case and empty values.

    Handles both the GFF3 form (key=value;) and the GTF2 form (key "value";).
    Returned in file order, so the attribute column can be rebuilt unchanged.
    """
    pairs = []
    for field in str(attributes).split(";"):
        field = field.strip()
        if not field:
            continue

        gtf2 = GTF2_PAIR.match(field)
        if gtf2:
            pairs.append((gtf2.group("key"), gtf2.group("value")))
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            pairs.append((key.strip(), value.strip()))
        else:
            # A bare field with no value at all; keep it so the round trip does
            # not silently drop information.
            pairs.append((field, ""))

    return pairs


def format_attributes(pairs):
    """Render (key, value) pairs back into a GFF3 attribute column."""
    return "".join("%s=%s;" % (key, value) for key, value in pairs)


def normalize_gff3_attribute_keys(pairs):
    """Return pairs whose keys follow GFF3's reserved-name convention.

    Values and pair order are deliberately untouched.  This is an explicit
    output-normalisation helper rather than part of :func:`split_attributes`,
    because readers must continue to preserve the spelling found in legacy
    files.
    """

    return [
        (GFF3_RESERVED_ATTRIBUTES.get(key.lower(), key.lower()), value)
        for key, value in pairs
    ]


def parse_attributes(attributes, lowercase_keys=True):
    """Attributes as a dict.

    Keys are lowercased by default, because the scripts look them up in lower
    case regardless of how the file spells them. Values are never altered. The
    first occurrence of a repeated key wins.
    """
    parsed = {}
    for key, value in split_attributes(attributes):
        parsed.setdefault(key.lower() if lowercase_keys else key, value)
    return parsed


def first_attribute(parsed, *keys, default=""):
    """The value of the first key present, in the order given."""
    for key in keys:
        if key in parsed:
            return parsed[key]
    return default


def require_attribute(parsed, key, context=""):
    """Fetch an attribute or stop with an explanation naming the record."""
    if key not in parsed:
        where = f" in {context}" if context else ""
        sys.exit(f"Missing '{key}' attribute{where}. Check your annotation.")
    return parsed[key]


def replace_attribute(pairs, key, value):
    """Set one attribute, matched case-insensitively, leaving the order intact.

    Returns a new list; the key is appended when it was not already present.
    """
    updated = []
    replaced = False
    for existing_key, existing_value in pairs:
        if existing_key.lower() == key.lower():
            updated.append((existing_key, value))
            replaced = True
        else:
            updated.append((existing_key, existing_value))

    if not replaced:
        updated.append((key, value))
    return updated

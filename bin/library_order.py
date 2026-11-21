"""Natural sort order for library names, shared by the report scripts.

Plain string sorting puts 100ng before 10ng and 1ng. natural_key compares the numbers
inside a name by value, so 1ng < 10ng < 100ng. A number followed by a mass unit (pg, ng,
ug, µg) is converted to ng first, so 100pg sorts before 1ng. Only a number with a unit
can have a decimal point, so 0.1ng is a tenth of a nanogram while S1.10 sorts after S1.2.
A unit must end the word, so 1ugly is the number 1 followed by text.
subworkflows/depth_dependent_analysis.nf applies the same rule to the Nextflow log.
"""
import re

# Either a mass (decimal allowed, unit not followed by a letter) or a plain integer. Both the
# micro sign (U+00B5) and the Greek mu (U+03BC) spell µg.
_NUMBER = re.compile(r'(\d+(?:\.\d+)?)\s*(pg|ng|ug|µg|μg)(?![a-z])|(\d+)', re.IGNORECASE)
_TO_NG = {'pg': 1e-3, 'ng': 1.0, 'ug': 1e3, 'µg': 1e3, 'μg': 1e3}


def natural_key(name):
    """Sort key: text pieces compare case-insensitively, numbers by (unit-scaled) value."""
    name = str(name)
    key, pos = [], 0
    for m in _NUMBER.finditer(name):
        if m.start() > pos:
            key.append((0, name[pos:m.start()].lower(), 0.0))
        if m.group(1):
            key.append((1, '', float(m.group(1)) * _TO_NG[m.group(2).lower()]))
        else:
            key.append((1, '', float(m.group(3))))
        pos = m.end()
    if pos < len(name):
        key.append((0, name[pos:].lower(), 0.0))
    return key


def natural_sorted(names):
    return sorted(names, key=natural_key)


def masses(name):
    """Every mass in a name as (value in ng, text as written), such as (0.1, '100pg')."""
    return [(float(m.group(1)) * _TO_NG[m.group(2).lower()], m.group(0))
            for m in _NUMBER.finditer(str(name)) if m.group(1)]

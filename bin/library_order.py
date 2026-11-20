"""Natural sort order for library names, shared by the report scripts.

Plain string sorting puts 100ng before 10ng and 1ng. natural_key compares the numbers
inside a name by value, so 0.1ng < 1ng < 10ng < 100ng. A number followed by a mass
unit (pg, ng, ug, µg) is converted to ng first, so 100pg sorts before 1ng.
lib/per_locus.nf applies the same rule to the Nextflow log.
"""
import re

_NUMBER = re.compile(r'(\d+(?:\.\d+)?)\s*(pg|ng|ug|µg)?', re.IGNORECASE)
_TO_NG = {'pg': 1e-3, 'ng': 1.0, 'ug': 1e3, 'µg': 1e3}


def natural_key(name):
    """Sort key: text pieces compare case-insensitively, numbers by (unit-scaled) value."""
    name = str(name)
    key, pos = [], 0
    for m in _NUMBER.finditer(name):
        if m.start() > pos:
            key.append((0, name[pos:m.start()].lower(), 0.0))
        unit = (m.group(2) or '').lower()
        key.append((1, '', float(m.group(1)) * _TO_NG.get(unit, 1.0)))
        pos = m.end()
    if pos < len(name):
        key.append((0, name[pos:].lower(), 0.0))
    return key


def natural_sorted(names):
    return sorted(names, key=natural_key)

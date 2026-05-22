"""Conversion type for all 12 possible single-base substitutions."""

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
_CONVERSIONS = {}
_ALL_TAGS = []

def _register(from_base, to_base):
    key = f"{from_base}>{to_base}"
    tag = f"{from_base.lower()}{to_base.lower()}"
    obj = {
        "tag": tag,
        "key": key,
        "from_base": from_base.upper(),
        "to_base": to_base.upper(),
        "comp_from": _COMPLEMENT[from_base.upper()],
        "comp_to": _COMPLEMENT[to_base.upper()],
    }
    _CONVERSIONS[key] = obj
    _CONVERSIONS[tag] = obj
    _ALL_TAGS.append(tag)
    return obj

# Register all 12 possible conversions
_register("C", "T")   # ct  - bisulfite
_register("T", "C")   # tc
_register("A", "G")   # ag
_register("G", "A")   # ga
_register("A", "C")   # ac
_register("C", "A")   # ca
_register("G", "T")   # gt
_register("T", "G")   # tg
_register("A", "T")   # at
_register("T", "A")   # ta
_register("C", "G")   # cg
_register("G", "C")   # gc


def get_conversion(value):
    """Resolve a conversion type from string key or tag.
    
    Accepts: 'C>T', 'c>T', 'C->T', 'c-to-t', 'ct', 'CT', etc.
    """
    if isinstance(value, dict):
        return value
    s = str(value).lower().replace(" ", "")
    # Normalize separators
    for sep in ["->", ">", "_", "-"]:
        if sep in s:
            parts = s.split(sep, 1)
            if len(parts) == 2:
                a, b = parts[0].strip()[:1], parts[1].strip()[:1]
                key = f"{a.upper()}>{b.upper()}"
                if key in _CONVERSIONS:
                    return _CONVERSIONS[key]
    # Try as a two-char tag
    if len(s) >= 2:
        tag = s[:2]
        if tag in _CONVERSIONS:
            return _CONVERSIONS[tag]
        # Try reversed char order
        tag_rev = s[1] + s[0]
        if tag_rev in _CONVERSIONS:
            return _CONVERSIONS[tag_rev]
    # Try 2-letter format (ct, ag, etc.)
    if s in _CONVERSIONS:
        return _CONVERSIONS[s]
    raise ValueError(
        f"Unknown conversion type: {value}. "
        f"Valid formats: 'C>T', 'ct', 'a2g', 'C->T'. "
        f"All 12 supported: {list(_CONVERSIONS.keys())}"
    )


ALL_CONVERSION_TAGS = list(_ALL_TAGS)


def list_conversions():
    """Return dict of all 12 conversion tags to descriptions."""
    desc = {
        "ct": "Bisulfite (C-to-T)",
        "tc": "T-to-C",
        "ag": "A-to-G",
        "ga": "G-to-A",
        "ac": "A-to-C",
        "ca": "C-to-A",
        "gt": "G-to-T",
        "tg": "T-to-G",
        "at": "A-to-T",
        "ta": "T-to-A",
        "cg": "C-to-G",
        "gc": "G-to-C",
    }
    return desc


class ConversionType:
    """Unified conversion type for any single-base substitution.
    
    Usage:
        conv = ConversionType("C>T")
        conv = ConversionType("ct")
        print(conv.target_base)      # C
        print(conv.converted_base)   # T
    """

    def __init__(self, value):
        self._data = get_conversion(value)

    @property
    def target_base(self):
        return self._data["from_base"]

    @property
    def converted_base(self):
        return self._data["to_base"]

    @property
    def complement_target(self):
        return self._data["comp_from"]

    @property
    def complement_converted(self):
        return self._data["comp_to"]

    @property
    def tag(self):
        return self._data["tag"]

    @property
    def key(self):
        return self._data["key"]

    def description(self):
        return f"{self.target_base}-to-{self.converted_base}"

    def __str__(self):
        return self.key

    def __repr__(self):
        return f"ConversionType('{self.key}')"

    def __eq__(self, other):
        if isinstance(other, ConversionType):
            return self.key == other.key
        return str(self) == str(other)

    def __hash__(self):
        return hash(self.key)

    @classmethod
    def from_str(cls, s):
        return cls(s)

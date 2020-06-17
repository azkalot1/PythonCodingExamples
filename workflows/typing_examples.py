from typing import List
from typing import Optional
import Path


def cook_foods(raw_foods: List[str]) -> List[str]:
    return [food.replace('raw', 'cooked') for food in raw_foods]


def str_or_none(opt_string: Optional[str] = None) -> Optional[str]:
    return opt_string


class TextFile:
    # Add type hints to TextFile"s __init__() method
    def __init__(self, name: str) -> None:
        self.text = Path(name).read_text()

    # Type annotate TextFile"s get_lines() method
    def get_lines(self) -> List[str]:
        return self.text.split("\n")


class MatchFinder:
    # Add type hints to __init__()'s strings argument
    def __init__(self, strings: List[str]) -> None:
        self.strings = strings

    # Type annotate get_matches()'s query argument
    def get_matches(self, query: Optional[str] = None) -> List[str]:
        return [s for s in self.strings if query in s] if query else self.strings

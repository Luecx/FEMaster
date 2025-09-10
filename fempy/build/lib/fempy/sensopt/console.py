# sensor_opt/console.py
from __future__ import annotations
from typing import Iterable

def hdr(title: str, width: int = 78) -> None:
    bar = "=" * width
    print("\n" + bar)
    print(title.strip())
    print(bar)

def subhdr(title: str, width: int = 78) -> None:
    bar = "-" * width
    print("\n" + title.strip())
    print(bar)

def kv(key: str, value, key_w: int = 24) -> None:
    k = (key + ":").ljust(key_w)
    print(f"  {k} {value}")

def bullet(items: Iterable[str], indent: int = 2, prefix: str = "• ") -> None:
    sp = " " * indent
    for s in items:
        print(f"{sp}{prefix}{s}")

def ok(msg: str) -> None:
    print(f"  [OK] {msg}")

def warn(msg: str) -> None:
    print(f"  [WARN] {msg}")

def err(msg: str) -> None:
    print(f"  [ERR] {msg}")

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import IO, Iterator, List, Optional, Type, Union

from more_itertools import peekable
from xopen import xopen

__all__ = [
    "ParsingError",
    "SAMHD",
    "SAMHeader",
    "SAMItem",
    "SAMReader",
    "SAMSQ",
    "read_sam",
]


class ParsingError(Exception):
    """
    Parsing error.
    """

    def __init__(self, line_number: int):
        """
        Parameters
        ----------
        line_number
            Line number.
        """
        super().__init__(f"Line number {line_number}.")
        self._line_number = line_number

    @property
    def line_number(self) -> int:
        """
        Line number.

        Returns
        -------
        Line number.
        """
        return self._line_number


@dataclasses.dataclass
class SAMHeader:
    """
    SAM header.

    Attributes
    ----------
    hd
        File-level metadata. Optional. If present, there must be only one
        @HD line and it must be the first line of the file.
    sq
        Reference sequence dictionary. The order of @SQ lines defines the
        alignment sorting order.
    rg
        Read group. Unordered multiple @RG lines are allowed.
    """

    hd: Optional[SAMHD] = None
    sq: List[SAMSQ] = dataclasses.field(default_factory=lambda: [])
    rg: List[str] = dataclasses.field(default_factory=lambda: [])


@dataclasses.dataclass
class SAMItem:
    """
    SAM item.

    Attributes
    ----------
    qname
        Query template NAME.
    """

    qname: str
    flag: str
    rname: str
    pos: str
    mapq: str
    cigar: str
    rnext: str
    pnext: str
    tlen: str
    seq: str
    qual: str
    remain: List[str]

    @classmethod
    def parse(cls: Type[SAMItem], line: str, line_number: int) -> SAMItem:
        values = line.strip().split("\t")
        try:
            args = tuple(values[:11]) + (values[11:],)
            item = cls(*args)
        except Exception:
            raise ParsingError(line_number)

        return item

    def copy(self) -> SAMItem:
        """
        Copy of itself.

        Returns
        -------
        SAM item.
        """
        from copy import copy

        return copy(self)


@dataclasses.dataclass
class SAMHD:
    vn: str
    so: Optional[str] = None

    @classmethod
    def parse(cls: Type[SAMHD], line: str, line_number: int) -> SAMHD:
        hd = cls(vn="")
        fields = line.strip().split("\t")

        try:
            assert fields[0] == "@HD"

            for f in fields[1:]:
                key, val = f.split(":")
                if key == "VN":
                    hd.vn = val
                elif key == "SO":
                    hd.so = val

            assert hd.vn != ""
        except Exception:
            raise ParsingError(line_number)

        return hd


@dataclasses.dataclass
class SAMSQ:
    sn: str
    ln: str

    @classmethod
    def parse(cls: Type[SAMSQ], line: str, line_number: int) -> SAMSQ:
        sq = cls("", "")
        fields = line.strip().split("\t")

        assert fields[0] == "@SQ"

        try:
            for f in fields[1:]:
                key, val = f.split(":")
                if key == "SN":
                    sq.sn = val
                elif key == "LN":
                    sq.ln = val

            assert sq.sn != ""
            assert sq.ln != ""
        except Exception:
            raise ParsingError(line_number)

        return sq


class SAMReader:
    """
    SAM reader.
    """

    def __init__(self, file: Union[str, Path, IO[str]]):
        """
        Parameters
        ----------
        file
            File path or IO stream.
        """
        if isinstance(file, str):
            file = Path(file)

        if isinstance(file, Path):
            file = xopen(file, "r")

        self._file = file
        self._lines = peekable(line for line in file)
        self._line_number = 0

        self._header = SAMHeader()
        try:
            next_line: str = self._lines.peek()
        except StopIteration:
            return

        while next_line.startswith("@"):

            line = self._next_line()

            if line.startswith("@HD"):
                self._header.hd = SAMHD.parse(line, self._line_number)
            elif line.startswith("@SQ"):
                self._header.sq.append(SAMSQ.parse(line, self._line_number))

            try:
                next_line = self._lines.peek()
            except StopIteration:
                break

    def read_item(self) -> SAMItem:
        """
        Get the next item.

        Returns
        -------
        Next item.
        """
        line = self._next_line()
        return SAMItem.parse(line, self._line_number)

    def read_items(self) -> List[SAMItem]:
        """
        Get the list of all items.

        Returns
        -------
        List of all items.
        """
        return list(self)

    def close(self):
        """
        Close the associated stream.
        """
        self._file.close()

    @property
    def header(self) -> SAMHeader:
        """
        File header.

        Returns
        -------
        Header.
        """
        return self._header

    def _next_line(self) -> str:
        line = next(self._lines)
        self._line_number += 1
        return line

    def __iter__(self) -> Iterator[SAMItem]:
        while True:
            try:
                yield self.read_item()
            except StopIteration:
                return

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        del exception_type
        del exception_value
        del traceback
        self.close()

    def __str__(self) -> str:
        return str(self._header)


def read_sam(file: Union[str, Path, IO[str]]) -> SAMReader:
    """
    Open a SAM file for reading.

    Parameters
    ----------
    file
        File path or IO stream.

    Returns
    -------
    SAM reader.
    """
    return SAMReader(file)

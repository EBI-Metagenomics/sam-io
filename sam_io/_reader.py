from __future__ import annotations

from typing import TypeVar, Type
import dataclasses
from pathlib import Path
from typing import IO, Iterator, List, Union, Optional, Dict

from more_itertools import peekable
from xopen import xopen

__all__ = ["ParsingError", "SAMItem", "SAMReader", "read_sam"]


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
    def parse(cls: Type[SAMHD], line: str) -> SAMHD:
        hd = cls(vn="")
        fields = line.strip().split("\t")

        assert fields[0] == "@HD"

        for f in fields[1:]:
            key, val = f.split(":")
            if key == "VN":
                hd.vn = val
            elif key == "SO":
                hd.so = val

        assert hd.vn != ""
        return hd


@dataclasses.dataclass
class SAMSQ:
    sn: str
    ln: str

    @classmethod
    def parse(cls: Type[SAMSQ], line: str) -> SAMSQ:
        sq = cls("", "")
        fields = line.strip().split("\t")

        assert fields[0] == "@SQ"

        for f in fields[1:]:
            key, val = f.split(":")
            if key == "SN":
                sq.sn = val
            elif key == "LN":
                sq.ln = val

        assert sq.sn != ""
        assert sq.ln != ""
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

        next_line: str = self._lines.peek()
        self._header = SAMHeader()
        while next_line.startswith("@"):

            line = self._next_line()

            if line.startswith("@HD\t"):
                self._header.hd = SAMHD.parse(line)
            elif line.startswith("@SQ\t"):
                self._header.sq.append(SAMSQ.parse(line))

            next_line = self._lines.peek()

    def read_item(self) -> SAMItem:
        """
        Get the next item.

        Returns
        -------
        Next item.
        """
        defline = self._next_defline()
        sequence = self._next_sequence()
        return SAMItem(defline, sequence)

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

    def _next_defline(self) -> str:
        while True:
            line = self._next_line()
            self._line_number += 1
            if line == "":
                raise StopIteration

            line = line.strip()
            if line.startswith(">"):
                return line[1:]
            if line != "":
                raise ParsingError(self._line_number)

    def _next_sequence(self) -> str:
        lines = []
        while True:
            line = self._next_line()
            if line == "":
                raise ParsingError(self._line_number)

            line = line.strip()
            if not line.startswith(">"):
                lines.append(line)
                if self._sequence_continues():
                    continue
                return "".join(lines)
            if line != "":
                raise ParsingError(self._line_number)

    def _sequence_continues(self):
        try:
            next_line = self._lines.peek()
        except StopIteration:
            return False

        if next_line == "":
            return False
        next_line = next_line.strip()
        return len(next_line) > 0 and not next_line.startswith(">")

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

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import IO, Iterator, List, Optional, Type, Union

from more_itertools import peekable
from xopen import xopen

__all__ = [
    "ParsingError",
    "SAMFlag",
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


class SAMFlag:
    """
    Bitwise flags.
    """

    def __init__(self, flag: int):
        self._flag = flag

    @property
    def value(self) -> int:
        return self._flag

    @property
    def read_paired(self) -> bool:
        return self._flag & 0x001 != 0

    @property
    def read_mapped_in_proper_pair(self) -> bool:
        return self._flag & 0x002 != 0

    @property
    def read_unmapped(self) -> bool:
        return self._flag & 0x004 != 0

    @property
    def mate_unmapped(self) -> bool:
        return self._flag & 0x008 != 0

    @property
    def read_reverse_strand(self) -> bool:
        return self._flag & 0x010 != 0

    @property
    def mate_reverse_strand(self) -> bool:
        return self._flag & 0x020 != 0

    @property
    def first_in_pair(self) -> bool:
        return self._flag & 0x040 != 0

    @property
    def second_in_pair(self) -> bool:
        return self._flag & 0x080 != 0

    @property
    def secondary_alignment(self) -> bool:
        return self._flag & 0x100 != 0

    @property
    def read_fails_filters(self) -> bool:
        return self._flag & 0x200 != 0

    @property
    def read_is_pcr_or_optical_duplicate(self) -> bool:
        return self._flag & 0x400 != 0

    @property
    def supplementary_alignment(self) -> bool:
        return self._flag & 0x800 != 0

    def __str__(self):
        return str(self._flag)

    def __repr__(self):
        return str(self._flag)


@dataclasses.dataclass
class SAMItem:
    """
    SAM item.

    Attributes
    ----------
    qname
        Query template NAME. Reads/segments having identical QNAME are regarded to come
        from the same template. A QNAME `*` indicates the information is unavailable.
    flag
        Combination of bitwise FLAGs.
    rname
        Reference sequence name of the alignment.
    pos
        1-based leftmost mapping position of the first CIGAR operation that "consumes" a
        reference base.
    mapq
        Mapping quality. It equals `−10 log10 Pr{mapping position is wrong}`, rounded to
        the nearest integer. A value 255 indicates that the mapping quality is not
        available.
    cigar
        CIGAR string.
    rnext
        Reference sequence name of the primary alignment of the next read in the
        template.
    pnext
        1-based position of the primary alignment of the next read in the template. Set
        as 0 when the information is unavailable.
    tlen
        Signed observed template length.
    seq
        Segment sequence. This field can be a `*` when the sequence is not stored.
    qual
        ASCII of base quality plus 33 (same as the quality string in the Sanger FASTQ
        format).
    remain
        Remaning fields not defined by SAM format.

    References
    ----------
    .. [SAMv1] https://samtools.github.io/hts-specs/SAMv1.pdf
    """

    qname: str
    flag: SAMFlag
    rname: str
    pos: int
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
            item = cls(
                values[0],
                SAMFlag(int(values[1])),
                values[2],
                int(values[3]),
                values[4],
                values[5],
                values[6],
                values[7],
                values[8],
                values[9],
                values[10],
                values[11:],
            )
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

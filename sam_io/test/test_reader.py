import os
import tarfile
from pathlib import Path

import importlib_resources as pkg_resources
import pytest

import sam_io


def test_read_example1(data_dir: Path):
    os.chdir(data_dir)

    with sam_io.read_sam("example1.sam") as file:
        header = file.header
        assert header.hd.vn == "1.6"
        assert header.hd.so == "queryname"
        assert len(header.sq) == 2
        assert header.sq[0].sn == "1"
        assert header.sq[0].ln == "1195517"
        assert header.sq[1].sn == "2"
        assert header.sq[1].ln == "6358"
        assert len(header.rg) == 0

        item = file.read_item()
        assert item.qname == "0b3cc96a-72f6-4d1e-a4ab-a43fa4993d34"
        assert item.flag == "0"
        assert item.rname == "2"
        assert item.pos == "1056"
        assert item.mapq == "60"
        assert item.cigar[:5] == "50M1D"
        assert item.cigar[-5:] == "1I38M"
        assert item.rnext == "*"
        assert item.pnext == "0"
        assert item.tlen == "0"
        assert item.seq[:5] == "TGTAG"
        assert item.seq[-5:] == "TCGAT"

        file.read_item()
        file.read_item()
        file.read_item()
        file.read_item()
        file.read_item()
        item = file.read_item()

        assert item.qname == "1b1b14f0-3bc2-465c-8259-053924e5cbb1"
        assert item.flag == "2064"
        assert item.seq[-5:] == "TCTGA"

        with pytest.raises(StopIteration):
            file.read_item()


def test_read_example2(data_dir: Path):
    os.chdir(data_dir)

    with sam_io.read_sam("example2.sam") as file:
        header = file.header
        assert header.hd is None
        assert len(header.sq) == 0
        assert len(header.rg) == 0

    items = list(sam_io.read_sam("example2.sam"))
    assert len(items) == 7

    assert items[0].qname == "0b3cc96a-72f6-4d1e-a4ab-a43fa4993d34"
    assert items[0].flag == "0"
    assert items[0].rname == "2"
    assert items[0].pos == "1056"
    assert items[0].mapq == "60"
    assert items[0].cigar[:5] == "50M1D"
    assert items[0].cigar[-5:] == "1I38M"
    assert items[0].rnext == "*"
    assert items[0].pnext == "0"
    assert items[0].tlen == "0"
    assert items[0].seq[:5] == "TGTAG"
    assert items[0].seq[-5:] == "TCGAT"

    assert items[6].qname == "1b1b14f0-3bc2-465c-8259-053924e5cbb1"
    assert items[6].flag == "2064"
    assert items[6].seq[-5:] == "TCTGA"


def test_read_example3(data_dir: Path):
    os.chdir(data_dir)

    with sam_io.read_sam("example3.sam") as file:
        header = file.header
        assert header.hd.vn == "1.6"
        assert header.hd.so == "queryname"
        assert len(header.sq) == 2
        assert header.sq[0].sn == "1"
        assert header.sq[0].ln == "1195517"
        assert header.sq[1].sn == "2"
        assert header.sq[1].ln == "6358"
        assert len(header.rg) == 0

        with pytest.raises(StopIteration):
            file.read_item()

    items = sam_io.read_sam("example3.sam").read_items()
    assert len(items) == 0


def test_read_example4(data_dir: Path):
    os.chdir(data_dir)

    with sam_io.read_sam("example4.sam") as file:
        header = file.header
        assert header.hd.vn == "1.6"
        assert header.hd.so == "queryname"
        assert len(header.sq) == 2
        assert header.sq[0].sn == "1"
        assert header.sq[0].ln == "1195517"
        assert header.sq[1].sn == "2"
        assert header.sq[1].ln == "6358"
        assert len(header.rg) == 0

        item = file.read_item()
        assert item.qname == "0b3cc96a-72f6-4d1e-a4ab-a43fa4993d34"
        assert item.flag == "0"
        assert item.rname == "2"
        assert item.pos == "1056"
        assert item.mapq == "60"
        assert item.cigar[:5] == "50M1D"
        assert item.cigar[-5:] == "1I38M"
        assert item.rnext == "*"
        assert item.pnext == "0"
        assert item.tlen == "0"
        assert item.seq[:5] == "TGTAG"
        assert item.seq[-5:] == "TCGAT"

        item = file.read_item()
        assert item.qname == "0bb88ec1-5961-48fa-b1b9-f5e72f116d68"

        file.read_item()
        file.read_item()
        file.read_item()
        file.read_item()
        file.read_item()

        with pytest.raises(StopIteration):
            file.read_item()

    items = sam_io.read_sam("example4.sam").read_items()
    assert len(items) == 7


def test_read_example5(data_dir: Path):
    os.chdir(data_dir)

    file = sam_io.read_sam("example5.sam")
    header = file.header
    assert header.hd is None
    assert len(header.sq) == 0
    assert len(header.rg) == 0

    with pytest.raises(StopIteration):
        file.read_item()
    file.close()

    items = sam_io.read_sam("example5.sam").read_items()
    assert len(items) == 0


def test_read_corrupted1(data_dir: Path):
    os.chdir(data_dir)

    with sam_io.read_sam("corrupted1.sam") as file:
        with pytest.raises(sam_io.ParsingError):
            file.read_item()


def test_read_corrupted2(data_dir: Path):
    os.chdir(data_dir)

    with pytest.raises(sam_io.ParsingError):
        sam_io.read_sam("corrupted2.sam")


def test_read_corrupted3(data_dir: Path):
    os.chdir(data_dir)

    file = sam_io.read_sam("corrupted3.sam")
    with pytest.raises(sam_io.ParsingError):
        file.read_item()


def test_read_corrupted4(data_dir: Path):
    os.chdir(data_dir)

    file = sam_io.read_sam("corrupted4.sam")
    file.read_item()
    with pytest.raises(sam_io.ParsingError):
        file.read_item()


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    os.chdir(tmp_path)

    content = pkg_resources.read_binary(sam_io.test, "examples.tar.gz")

    with open("examples.tar.gz", "wb") as f:
        f.write(content)

    tar = tarfile.open("examples.tar.gz", "r:gz")
    tar.extractall()
    tar.close()

    os.unlink("examples.tar.gz")

    return tmp_path

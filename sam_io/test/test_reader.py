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

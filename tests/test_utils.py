import pytest
from dataclasses import dataclass

from betterbasequals import utils

@dataclass
class DummyPilupColumn:
    """Class for testing pileup."""
    reference_name: str
    reference_pos: int

def test_zip_pileups():
    IL1 = [1,2,4,5,7,8,9,11,13,15]
    cL1 = ['chr1','chr2','chr3','chr4']
    L1 = []
    tup1 = set()
    for c in cL1:
        for x in IL1:
            L1.append(DummyPilupColumn(c,x))
            tup1.add((c,x))
    IL2 = [1,2,4,5,7,9,10,12,13,14,15]
    cL2 = ['chr1','chr3','chr4']
    L2 = []
    tup2 = set()
    for c in cL2:
        for x in IL2:
            L2.append(DummyPilupColumn(c,x))
            tup2.add((c,x))
    IL3 = [1,3,5,7,8,9,10,13,14,15]
    cL3 = ['chr1','chr2','chr4']
    L3 = []
    tup3 = set()
    for c in cL3:
        for x in IL3:
            L3.append(DummyPilupColumn(c,x))
            tup3.add((c,x))
    G1 = (x for x in L1)
    G2 = (x for x in L2)
    G3 = (x for x in L3)
    common = tup1 & tup2 & tup3
    found_common = set()
    for x,y,z in utils.zip_pileups(G1,G2,G3):
        assert x.reference_pos == y.reference_pos == z.reference_pos and x.reference_name == y.reference_name == z.reference_name
        found_common.add((x.reference_name, x.reference_pos))
    assert common == found_common

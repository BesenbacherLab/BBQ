import csv
from itertools import product, count
from operator import add
from functools import reduce, partial

complement = str.maketrans("ATCGN", "TAGCN")


def reverse_complement(s):
    return s.translate(complement)[::-1]


def regions(bed_file):
    with open(bed_file) as fp:
        reader = csv.reader(fp, delimiter="\t")
        for line in reader:
            yield line[0], int(line[1]), int(line[2])

def regions_if_sorted(bed_file):
    with open(bed_file) as fp:
        reader = csv.reader(fp, delimiter="\t")
        for line in reader:
            yield line[0], int(line[1]), int(line[2])

def generate_mutation_types(k):
    if k % 2 == 0:
        raise ValueError("k must be uneven")
    types = list()
    mut_pos = (k - 1) // 2
    for mer in map(
        partial(reduce, add),
        product(*(["ATGC"] * mut_pos + ["TC"] + ["ATGC"] * mut_pos)),
    ):
        alts = "ATG"
        if mer[mut_pos] == "T":
            alts = "AGC"
        for alt in alts:
            after = mer[:mut_pos] + alt + mer[mut_pos + 1 :]
            mut = ">".join([mer[mut_pos],alt])
            types.append("_".join([mut, mer]))
    return sorted(types)

mutation_types_3_mer = dict(zip(generate_mutation_types(3), count()))

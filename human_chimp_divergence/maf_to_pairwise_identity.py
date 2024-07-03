# Code adapted from https://github.com/makovalab-psu/great-ape-Y-evolution
"""
Read alignments in maf format and output pairwise identity stats.
"""

from gzip import open as gzip_open
import copy
from io import StringIO

from .maf_reader import maf_alignments
from .fasta_file import FastaFile

def maf_to_pairwise_identity(fasta_files, input_maf, identity_file, chimp_vcf=None):
    reportProgress = 1

    pairToData = {}  # base-by-base alignment events
    nameToLengths = {}  # species.contig lengths as reported in maf file

    mafBlockNumber = 0

    matchEvent = 4  # match
    subEvent = 3  # substitution
    insEvent = 2  # inserted in ref
    delEvent = 1  # deleted from ref

    eventToCh = {matchEvent: "M", subEvent: "S", insEvent: "I", delEvent: "D"}

    # open fasta files
    # nameToLookup = {}
    # for name in fasta_files:
    #     filename = fasta_files[name]

    nameToLookup = {}
    filename = fasta_files['hg19']
    if (filename.endswith(".gz")) or (filename.endswith(".gzip")):
        nameToLookup['hg19'] = FastaFile(gzip_open(filename, "rt").readlines())
    else:
        nameToLookup['hg19'] = FastaFile(open(filename, 'r').readlines())

    with open(input_maf, "r") as file:
        maf_data = file.read()
    maf_data2 = StringIO(maf_data)

    for a in maf_alignments(maf_data2, yieldHeader=True, discardWeeds=False):
        if isinstance(a, str):  # this is a maf header
            continue

        mafBlockNumber += 1
        # if reportProgress != None:
        #     if (mafBlockNumber == 1) or (mafBlockNumber % reportProgress == 0):
        #         print("progress: maf block %d (line %d)" % (mafBlockNumber, a.lineNum), file=sys.stderr)

        for c1, c2 in pairwise_alignments(a.block):
            if c1.ref not in nameToLengths:
                nameToLengths[c1.ref] = {}
            if c1.contig not in nameToLengths[c1.ref]:
                nameToLengths[c1.ref][c1.contig] = c1.srcSize
            if c2.ref not in nameToLengths:
                nameToLengths[c2.ref] = {}
            if c2.contig not in nameToLengths[c2.ref]:
                nameToLengths[c2.ref][c2.contig] = c2.srcSize

            accumulate_base_by_base(c1, c2, pairToData)

    # if reportProgress != None:
    #     print("progress: %d maf blocks total" % (mafBlockNumber), file=sys.stderr)

    # report stats
    refToNonNCount = {}  # counts non-Ns in a species

    # print("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
    #     % ("ref", "other", "refNonN", "m", "mm", "i", "d", "aligned", "identity"))

    refs = list(fasta_files.keys())

    for refName, otherName in [(refs[0], refs[1])]:
        speciesPair = (refs[0], refs[1])
        if speciesPair not in pairToData:
            print("%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA" % (refName, otherName))
            continue

        if refName not in refToNonNCount:
            refLookup = nameToLookup[refName]
            refContigSeen = set()
            refToNonNCount[refName] = 0
            for refContig in refLookup.keys():
                refContigSeen.add(refContig)
                refSeq = refLookup[refContig].upper()
                refToNonNCount[refName] += len([nuc for nuc in refSeq if (nuc != "N")])
                if (refContig in nameToLengths[refName]) and (len(refSeq) != nameToLengths[refName][refContig]):
                    assert False, ("maf input says len(%s.%s) is %d, but in fasta file it is %d"
                        % ( refName, refContig, nameToLengths[refName][refContig], refSeq))
            for refContig in nameToLengths[refName]:
                assert refContig in refContigSeen, (
                    "maf input refers to %s.%s, but fasta file has no such sequence"
                    % (refName, refContig)
                )

        nonNCount = refToNonNCount[refName]

        refData = pairToData[speciesPair]
        eventToCount = {
            matchEvent: 0,
            subEvent: 0,
            insEvent: 0,
            delEvent: 0,
        }
        for refContig in refData:
            contigEventVector = refData[refContig]
            refPositions = contigEventVector.keys()
            refPositions = sorted(refPositions)
            for refPos in refPositions:
                eventToCount[contigEventVector[refPos][0]] += 1

        m = eventToCount[matchEvent]
        mm = eventToCount[subEvent]
        i = eventToCount[insEvent]
        d = eventToCount[delEvent]

        # print( "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%"
        #     % (refName, otherName, nonNCount, m, mm, i, d,
        #         100.0 * (m + mm + i) / nonNCount,
        #         100.0 * m / (m + mm + i) ) )

        for refContig in refData:
            contigEventVector = refData[refContig]
            refPositions = contigEventVector.keys()
            refPositions = sorted(refPositions)
            with open(identity_file, "wt") as file:
                for refPos in refPositions:
                    file.write(
                        "%s\t%d\t%s\t%s\t%s\n"
                        % (
                            refContig,
                            refPos,
                            eventToCh[contigEventVector[refPos][0]],
                            contigEventVector[refPos][1],
                            contigEventVector[refPos][2],
                        )
                    )

            if chimp_vcf:
                with open(chimp_vcf, "wt") as file:
                    file.write(
                        '##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchimp\n'
                    )
                    for refPos in refPositions:
                        if (contigEventVector[refPos][2] == "N") | (contigEventVector[refPos][2] == "-"):
                            alt = "."
                            gt = "."
                        elif (contigEventVector[refPos][2] == contigEventVector[refPos][1]):
                            alt = "."
                            gt = "0"
                        else:
                            alt = contigEventVector[refPos][2]
                            gt = "1"

                        file.write(
                            "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT\t{}\n".format(
                                refContig, refPos + 1, contigEventVector[refPos][1], alt, gt))


# pairwise_alignments--
# yield reference-centric pairwise blocks from a maf alignment block
#
# The first component in each reported pairwise block is always on the "+"
# strand. Note that each block is reported twice, once with one species as the
# first component, once with the other as first.


def pairwise_alignments(block):
    if len(block) < 2:
        return

    # identify all the relevant components

    nameToComponents = {}
    speciesPresent = []

    for componentIx, c in enumerate(block):
        if c.ref not in nameToComponents:
            nameToComponents[c.ref] = [componentIx]
            speciesPresent += [c.ref]
        else:
            nameToComponents[c.ref] += [componentIx]

    if len(speciesPresent) < 2:
        return

    # yield the pairwise alignments

    componentToReverse = {}

    for refName, otherName in all_pairs(speciesPresent, ordered=True):
        for componentIx1 in nameToComponents[refName]:
            c1 = block[componentIx1]
            flipForStrand = False
            if c1.strand == "-":
                if componentIx1 not in componentToReverse:
                    componentToReverse[componentIx1] = reverse_of_component(c1)
                c1 = componentToReverse[componentIx1]
                flipForStrand = True

            for componentIx2 in nameToComponents[otherName]:
                c2 = block[componentIx2]
                if flipForStrand:
                    if componentIx2 not in componentToReverse:
                        componentToReverse[componentIx2] = reverse_of_component(c2)
                    c2 = componentToReverse[componentIx2]

                yield (c1, c2)


# reverse_of_component--
# reverse complement a component


def reverse_of_component(c):
    c = copy.copy(c)
    c.start = c.srcSize - (c.start + c.length)
    c.strand = "+" if (c.strand == "-") else "-"
    c.nucs = reverse_complement(c.nucs)
    return c


# accumulate_base_by_base--
# accumulate each reference base into a base-by-base vector of alignment
# events


def accumulate_base_by_base(cRef, cOther, pairToData):
    matchEvent = 4  # match
    subEvent = 3  # substitution
    insEvent = 2  # inserted in ref
    delEvent = 1  # deleted from ref

    assert cRef.strand == "+"

    refName = cRef.ref
    refContig = cRef.contig  # note that this can be None
    otherName = cOther.ref
    otherContig = cOther.contig  # note that this can be None

    refNucs = cRef.nucs.upper()
    otherNucs = cOther.nucs.upper()

    speciesPair = (refName, otherName)
    if speciesPair not in pairToData:
        pairToData[speciesPair] = {}
    refData = pairToData[speciesPair]

    if refContig not in refData:
        refData[refContig] = {}
    contigEventVector = refData[refContig]

    refPos = cRef.start

    n_ins = 0
    n_del = 0
    n_subst = 0
    n_match = 0

    for refCh, otherCh in zip(refNucs, otherNucs):
        refIsGap = refCh == "-"
        otherIsGap = otherCh == "-"
        if (refIsGap) and (otherIsGap):
            continue
        if refCh == "N":
            continue
        if refCh == otherCh:
            event = matchEvent
            n_match += 1
        elif refIsGap == otherIsGap:
            event = subEvent
            n_subst += 1
        elif otherIsGap:
            event = insEvent
            n_ins += 1
        else:
            event = delEvent
            n_del += 1

        temp = [event, refCh, otherCh]

        if refPos not in contigEventVector:
            contigEventVector[refPos] = temp
        elif event > contigEventVector[refPos][0]:
            contigEventVector[refPos] = temp

        if refCh != "-":
            refPos += 1


# all pairs--
# yield all pairs from a set
#
# if ordered is false, pairs (a,b) and (b,a) are considered the same, and only
# one is generated


def all_pairs(items, ordered=False):
    for ix in range(len(items) - 1):
        for iy in range(ix + 1, len(items)):
            yield (items[ix], items[iy])
            if ordered:
                yield (items[iy], items[ix])


# reverse_complement--


def reverse_complement(nucs):
    return nucs[::-1].translate(str.maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn", "TGCASWYRKMVHDBNtgcaswyrkmvhdbn"))

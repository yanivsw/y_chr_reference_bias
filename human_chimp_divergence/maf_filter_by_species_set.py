# Code adapted from https://github.com/makovalab-psu/great-ape-Y-evolution

"""
Read alignments in maf format and output those blocks that have a specified
set of species, and no other species.
"""

from .maf_reader import maf_alignments
from io import StringIO
import gzip

def filter_maf(human, chimp, input_maf, contigs, min_alignment_score, filtered_maf):
    discardDuplicates = False
    reportProgress    = 1

    nothingHasPrinted = True

    mafBlockNumber = 0
    mafBlocksKept = mafBlocksDiscarded = 0

    if (input_maf.endswith(".gz")) or (input_maf.endswith(".gzip")):
        with gzip.open(input_maf,'r') as file:        
            maf_data = file.read().decode("utf-8")
    else:
        with open(input_maf, 'r') as file:
            maf_data = file.read()

    maf_data2 = StringIO(maf_data)

    out_str = ''
    for a in maf_alignments(maf_data2,saveLines=True):
        mafBlockNumber += 1
        # if (reportProgress != None) and (mafBlockNumber % reportProgress == 0):
        # 	print("progress: maf block %d (line %d), %d kept, %d discarded" % (mafBlockNumber,a.lineNum,
        # 	                 mafBlocksKept,mafBlocksDiscarded), file=sys.stderr)

        score = float(a.score)
        if contigs['chimp'] is not None:
            keep_block = ( ( ((a.block[0].contig == contigs['human']) & (a.block[0].ref == human)) | ((a.block[1].contig == contigs['human']) & (a.block[1].ref == human) ) ) &
                           ( ((a.block[0].contig == contigs['chimp']) & (a.block[0].ref == chimp)) | ((a.block[1].contig == contigs['chimp']) & (a.block[1].ref == chimp) ) ) )
        else:
            keep_block = ( ( (a.block[0].contig == contigs['human']) & (a.block[0].ref == human) ) | ((a.block[1].contig == contigs['human']) & (a.block[1].ref == human) ) )
    
        keep_block = keep_block & (score >= min_alignment_score)

        if not keep_block:
            mafBlocksDiscarded += 1
            continue

        if (discardDuplicates):
            refs = set()
            hasDuplicate = False
            for c in a.block:
                if (c.ref in refs):
                    hasDuplicate = True
                    break
                refs.add(c.ref)
            if (hasDuplicate):
                mafBlocksDiscarded += 1
                continue

        mafBlocksKept += 1
        if (nothingHasPrinted):
            nothingHasPrinted = False
        else:
            out_str += ("\n")
        
        out_str += ("a score=" + a.score + "\n")
        for c in a.block:
            out_str += (c.line + '\n')

    with open(filtered_maf, 'w') as file:
        file.write(out_str)

    # if (reportProgress is not None):
    #     print("progress: %d maf blocks total" % (mafBlockNumber), file=sys.stderr)

    # if (mafBlocksKept + mafBlocksDiscarded == 0):
    #     print("no maf blocks kept, none discarded", file=sys.stderr)
    # else:
    #     print("%d maf blocks kept (%.2f%%), %d discarded" % (mafBlocksKept,
    #                      100.0*mafBlocksKept/(mafBlocksKept+mafBlocksDiscarded),
    #                      mafBlocksDiscarded), file=sys.stderr)

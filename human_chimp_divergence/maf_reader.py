# Code adapted from https://github.com/makovalab-psu/great-ape-Y-evolution

from sys import stderr

# maf_alignments--
#	yields the next maf alignment from a file
#
# The input is a maf file, as per the spec at
#   http://genome.ucsc.edu/FAQ/FAQformat.html#format5
#
# discardWeeds means that any component 'off in the weeds' is discarded. These
# are components that are off the end of the sequence. Normally these would
# cause an exception. This option was added because halToMaf (part of the
# progressiveCactus package) occasionally creates such components.

class Alignment:
    pass

class Component:
    pass

def maf_alignments(f,yieldHeader=False,saveLines=False,discardWeeds=False):
    mafHeader = None
    a = None

    context = []
    lineNum = 0
    discardedComponentCount = 0
    for line in f:
        lineNum += 1
        line = line.strip()

        if (line == "##eof maf"):
            continue
        if line.startswith("##"):
            if mafHeader is None:
                mafHeader = [line]
            elif isinstance(mafHeader, list):
                mafHeader.append(line)
            continue

        if mafHeader is None:
            mafHeader = ""
        elif isinstance(mafHeader, list):
            mafHeader = "\n".join(mafHeader)
            if yieldHeader:
                yield mafHeader

        if line.startswith("#"):
            #print line
            continue

        if (line == "") and (a is not None):
            if (discardedComponentCount != 0) and (a.block == []):
                a = None
                continue
            assert (a.block != []), \
                   "lacking components for alignment at line %d" \
                 % (a.lineNum)
            yield a
            a = None
            continue

        if line.startswith("a "):
            assert a is None, \
                   "unexpected score line (line %d)\n%s" \
                 % (lineNum,"\n".join(context))
            a = Alignment()
            a.score = line.split()[1].split("=")[1]
            a.block = []
            a.lineNum = lineNum
            discardedComponentCount = 0
            continue

        if line.startswith("a"):
            assert a is None, \
                   "unexpected \"a\" line (line %d)\n%s" \
                 % (lineNum,"\n".join(context))
            a = Alignment()
            a.score = None
            a.block = []
            a.lineNum = lineNum
            discardedComponentCount = 0
            continue

        if (not line.startswith("s ")):
            continue

        assert (a is not None), \
               "lacking score line before alignment line (line %d)\n%s" \
             % (lineNum,"\n".join(context))

        c = Component()

        try:
            fields = line.split()
            if len(fields) != 7:
                raise ValueError
            c.lineNum = lineNum
            if saveLines:
                c.line = line
            src = fields[1]
            if "." in src:
                c.ref, c.contig = src.split(".", 1)
            else:
                c.ref, c.contig = src, None
            c.start   = int(fields[2])
            c.length  = int(fields[3])
            c.strand  =     fields[4]
            c.srcSize = int(fields[5])
            c.nucs    =     fields[6]
            if (c.start  < 0):
                raise ValueError
            if (c.length < 0):
                raise ValueError
            if (c.strand not in ["+","-"]):
                raise ValueError
            if (c.srcSize < c.start+c.length):
                if (discardWeeds):
                    c = None
                    if (len(line) > 80):
                        line = line[:77]+"..."
                    print("WARNING: discarding line (%d); it is in the weeds (component %d):\n\"%s\"" \
                        % (lineNum,1+len(a.block)+discardedComponentCount,line), file=stderr)
                else:
                    raise ValueError
        except ValueError:
            assert (False), \
                   "bad line (%d):\n\"%s\"" % (lineNum,line)
            #assert (False), \
            #       "bad line (%d): %s %s %s %s\"%s\"" \
            #     % (lineNum,
            #        (c.start  < 0),
            #        (c.length < 0),
            #        (c.srcSize < c.start+c.length),
            #        (c.strand not in ["+","-"]),
            #        line)

        if (c is None):
            discardedComponentCount += 1
            continue

        if (a.block == []):
            alignmentColumns = len(c.nucs)
        else:
            assert (len(c.nucs) == alignmentColumns), \
                   "inconsistent alignment length, expected %d but got %d (line %d): \"%s\"" \
                 % (alignmentColumns,len(c.nucs),lineNum,line)

        a.block += [c]

    if (a is not None):
        if (discardedComponentCount != 0) and (a.block == []):
            pass
        else:
            assert (a.block != []), \
                   "lacking components for alignment at line %d" \
                 % (a.lineNum)
            yield a


# Code adapted from https://github.com/makovalab-psu/great-ape-Y-evolution

from struct   import *

from collections import MutableMapping

class FastaFile(MutableMapping):
    def __init__(self,f,unmask=False):
        self.debug = []

        # if isinstance(f, str):
        #     with open(f, 'rb') as file:
        #         f = file.readlines()

        self.file   = f
        self.unmask = unmask

        self.seqCount = 0
        self.index = {}
        self.indexByAlias = {}
        for (name,seq) in self.read_sequences(f):
            if (unmask):
                seq = seq.upper()
            self.index[name] = seq
            self.indexByAlias[self.seqCount] = seq
            self.seqCount += 1
        self.preloadedCount = self.seqCount
        if (self.seqCount == 1):
            self.indexByAlias[None] = self.indexByAlias[0]

    def __getitem__(self,name):
        if (name in self.index):
            return self.index[name]
        else:
            return self.indexByAlias[name]

    def __setitem__(self,name):
        return

    def __delitem__(self,name):
        return

    def __len__(self):
        return len(self.mylist)

    def __iter__(self):
        for i in self.mylist:
            yield i
    
    def keys(self):
        return self.index.keys()

    def read_sequences(self,f):
        seqName = None
        lineNum = 0
        for line in f:
            lineNum += 1
            line = line.strip()

            if (line.startswith(">")):
                if (seqName is not None):
                    seq = "".join(seq)
                    yield (seqName,seq)
                fields = line[1:].split()
                if (fields == []):
                    assert (False), \
                           "sequence has no name (at line %d)" % lineNum
                seqName = fields[0]
                seq = []
            elif (seqName is None):
                assert (False), "first sequence has no header"
            else:
                seq += [line]

        if (seqName is not None):
            seq = "".join(seq)
            yield (seqName,seq)

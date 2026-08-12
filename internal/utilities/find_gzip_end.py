#!/usr/bin/env python3
import sys
import zlib

path = sys.argv[1]
decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
consumed = 0
with open(path, "rb") as handle:
    while block := handle.read(4 * 1024 * 1024):
        consumed += len(block)
        decoder.decompress(block)
        if decoder.unused_data:
            print(consumed - len(decoder.unused_data))
            break
    else:
        print(consumed)

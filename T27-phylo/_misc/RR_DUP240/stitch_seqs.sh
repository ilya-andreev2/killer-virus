#!/bin/sh

for filename in DUP*.fasta; do
    echo $filename
done

for f in DUP*.fasta; do (cat "${f}"; echo) >> RR_DUP240.fsa; done
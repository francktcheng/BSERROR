#!/bin/bash

PATH1='*.pdf'
# PATH2='S_*.pdf'

for f in $PATH1
do
  pdfcrop --margins=5 $f $f
done


#!/bin/bash -ue
# visualise gfa - need to separate, this takes a while
vg view -F -p -d pggb.gfa  | dot -Tpng -o out.vg.png

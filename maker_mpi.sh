#!/bin/bash

mpiexec -n $NSLOTS /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/mpi_maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl

cd bccFe

## make perfect slab
## the command line inputs to make-slab are:
## [cell file] [infile] [Rcut]

../../make-slab cell_bccFe infile_bccedge_t0-11 300 > perf

## self-consistently evaluate edge dislocation geometry
## the command line inputs to anisotropic-xyz-ref are:
## [cell file] [infile] [undislocated xyz file] [reference file =  dislocated xyz file from previous iteration]

../../anisotropic-xyz-ref cell_bccFe infile_bccedge_t0-11 perf perf > edge_1
for i in $(seq 1 9); do
	../../anisotropic-xyz-ref cell_bccFe infile_bccedge_t0-11 perf edge_$(echo "$i") > edge_$(echo "$i+1" | bc)
done

## compute and output strain
#../../anisotropic-xyz-ref-outputstrain cell_bccFe infile_bccedge_t0-11 perf edge_10 strain > edge_11

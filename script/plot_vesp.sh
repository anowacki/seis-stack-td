#!/bin/bash
# Plot a vespagram.

usage() {
	cat <<-END >&2
	Usage: $(basename $0) (stack_vespa options)
	Plot a vespagram using stack_vespa.  Command line options are
	passed to the program, but do not use -o, as this is reserved
	for the plotting script.
	
	Usage for stack_vespa
	---------------------
	$(stack_vespa 2>&1)
	END
	exit 1
}

# Check we haven't tried to use the -o option
for a in "$@"; do
	[ "$a" = "-o" ] && usage
done

# Make temporary files
CPT=$(mktemp /tmp/plot_vesp.sh.cptXXXXXX)
FIG=$(mktemp /tmp/plot_vesp.sh.psXXXXXX)
GRD=$(mktemp /tmp/plot_vesp.sh.grdXXXXXX)
trap 'rm -f "$CPT" "$FIG" "$GRD"' EXIT

# Create vespagram, passing options.  Anything read from stdin will be passed in
stack_vespa -o "$GRD" "$@" 2>/dev/null || usage

# Get info from grid file
read t1 t2 dt s1 s2 ds <<< $(grdinfo "$GRD" |
	awk '/x_min/{print $3,$5,$7} /y_min/{print $3,$5,$7}')
read lt ls <<< $(echo $t1 $t2 $s2 $s2 | awk '{printf("%f %f", ($2-$1)/5, ($4-$3)/4)}')
makecpt -Z -Cpolar -T-1/1/0.1 > "$CPT"

grdimage "$GRD" -JX8c/6c -R$t1/$t2/$s1/$s2 -C"$CPT" -P \
	-Ba$lt":Time / s:"/a0.5":Slowness / s/deg:"":.Slowness vespagram:"nSeW > "$FIG" &&
gv "$FIG"

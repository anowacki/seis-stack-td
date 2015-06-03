#!/bin/bash
# Plot FK analysis.

usage() {
	cat <<-END >&2
	Usage: $(basename $0) (stack_fk options)
	Plot an array response using stack_fk.  Command line options are
	passed to the program, but do not use -o, as this is reserved
	for the plotting script.
	
	$(stack_fk)
	END
	exit 1
}

# Check we haven't tried to use the -o option
for a in "$@"; do
	[ "$a" = "-o" ] && usage
done

# Make temporary files
CPT=$(mktemp /tmp/plot_fk.sh.cptXXXXXX)
FIG=$(mktemp /tmp/plot_fk.sh.psXXXXXX)
GRD=$(mktemp /tmp/plot_fk.sh.grdXXXXXX)
trap 'rm -f "$CPT" "$FIG" "$GRD"' EXIT

# Plotting defaults
ls=2

# Create vespagram, passing options.  Anything read from stdin will be passed in
stack_fk -o "$GRD" "$@" 2>/dev/null || usage

# Get info from grid file
read smin smax ds <<< $(grdinfo "$GRD" | awk '/x_min/{print $3,$5,$7}')
makecpt -Z -Chaxby -I -T0/1/0.05 > "$CPT" 2>/dev/null

# Plot power
grdimage "$GRD" -JX8c/8c -R$smin/$smax/$smin/$smax -C"$CPT" -P -K \
	-Ba$ls":@%2%u@-x@-@%% / s/deg:"/a$ls":@%2%u@-y@-@%% / s/deg:"":.Beam power:"nSeW > "$FIG" &&
# Plot circles for slowness and lines for azimuth
awk -v smax=$smax 'BEGIN {
	pi = 4*atan2(1,1)
	for (r=2; r<=smax*sqrt(2); r+=2) {
		print ">"
		for (theta=0; theta<=2*pi; theta+=pi/180) print r*sin(theta), r*cos(theta)
	}
	for (theta=0; theta<2*pi; theta+=pi/6) {
		print ">"
		for (r=0; r<=100; r+=100) print r*sin(theta), r*cos(theta)
	}
	}' | psxy -J -R -m -W0.5p,- -O >> "$FIG"
gv "$FIG"

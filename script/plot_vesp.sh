#!/bin/bash
# Plot a vespagram.

usage() {
	cat <<-END >&2
	Usage: $(basename $0) (stack_vespa options)
	Plot a vespagram using stack_vespa.  Command line options are
	passed to the program, but do not use -o, as this is reserved
	for the plotting script.
	
	Options:
	   -annot [name] [time] [slow] : Add an annotation at time time and
	              slowness slow.  _s in name will be replaced by spaces
	              in plotting.
	   -phase [phase1(,phase2)]    : Plot phases using the event-array geography.
	              Quote the list and add "-mod [model]" to change model.
	
	Usage for stack_vespa
	---------------------
	$(stack_vespa 2>&1)
	END
	exit 1
}

die() { echo "$@"; exit 1; }

plot_cross() {
	# Plot a cross at a given time, slowness location
	# Usage: plot_cross time slowness
	[ $# -ne 2 ] && { echo "Error with plot_cross">&2; return 1; }
	echo "$1" "$2" | psxy -J -R -S+0.5c -W2.5p,white -O -K >> "$FIG"
	echo "$1" "$2" | psxy -J -R -S+0.45c -W1p,black -O -K >> "$FIG"
}

plot_annot_cross() {
	# Plot an annotated cross at a given time, slowness location
	# Usage: plot_annot_cross time slowness "label"
	[ $# -ne 3 ] && { echo "Error with plot_annot_cross">&2; return 1; }
	plot_cross "$1" "$2"
	echo "$1" "$2" 10 0 0 BL "$3" | pstext -J -R -D0.2c/0.2c -O -K >> "$FIG"
}

# Process arguments
while [ "$1" ]; do
	case "$1" in
		# Arguments the script knows about
		-phase) phases="$2"; shift 2;;
		-annot) name_list=("${name_list[@]}" "$2")
		        time_list=("${time_list[@]}" "$3")
		        slow_list=("${slow_list[@]}" "$4"); shift 4;;
		# Arguments passed to stack_fk are anything we don't know about
		*) break;;
	esac
done

# Check we're not trying to use to -o option
for arg in "$@"; do
	[ "$arg" = "-o" ] && { echo "Do not use option '-o' with plotting script"; usage; }
done
[ $# -eq 0 ] && usage

# Make temporary files
CPT=$(mktemp /tmp/plot_vesp.sh.cptXXXXXX)
FIG=$(mktemp /tmp/plot_vesp.sh.psXXXXXX)
GRD=$(mktemp /tmp/plot_vesp.sh.grdXXXXXX)
trap 'rm -f "$CPT" "$FIG" "$GRD"' EXIT

# Create vespagram, passing options.  Anything read from stdin will be passed in
stack_vespa -o "$GRD" "$@" || { echo "Error running stack_vespa" >&2; exit 1; }

# Get info from grid file
read gcarc baz evdp t1 t2 dt s1 s2 ds <<< $(grdinfo "$GRD" |
	awk '/Command:/{print $(NF-4), $(NF-2), $NF}
		/x_min/{print $3,$5,$7} /y_min/{print $3,$5,$7}')
read lt ls <<< $(echo $t1 $t2 $s2 $s2 | awk '{printf("%f %f", ($2-$1)/5, ($4-$3)/4)}')
makecpt -Z -Cpolar -T-1/1/0.1 > "$CPT"

# Plot grid file
grdimage "$GRD" -JX8c/6c -R$t1/$t2/$s1/$s2 -C"$CPT" -P \
	-Ba$lt":Time / s:"/a0.5":Slowness / s/deg:"":.Slowness vespagram:"nSeW -K > "$FIG" ||
	die "Error plotting grid file"

# Add phase arrivals using taup if available
if [ "$phases" ]; then
	command -v taup_time >/dev/null 2>&1 ||
		{ echo "Cannot find taup_time; no phases will be plotted"; break; }
	list=$(taup_time -ph $phases -h $evdp -deg $gcarc | awk 'NR>=6')
	echo "$list" | while read gcarc_taup evdp_taup phase time slowness takeoff incident \
			distance blank pure_name; do
		plot_annot_cross $time $slowness $phase
	done || die "Error plotting phase arrivals"
fi

# Add annotations if any
for ((i=0; i<${#name_list[@]}; i++)); do
	plot_annot_cross "${time_list[i]}" "${slow_list[i]}" "${name_list[i]}"
done || die "Error plotting annotations"

# Add plot information at the top
printf "%f %f 10 0 0 BL @~D@~ = %0.1f@~\\260@~  baz = %0.1f@~\\260@~\n" $t1 $s2 $gcarc $baz |
	pstext -J -R -D0c/0.2c -N -O -K >> "$FIG"

# Finish postscript
psxy -J -R -T -O >> "$FIG"

gv "$FIG"

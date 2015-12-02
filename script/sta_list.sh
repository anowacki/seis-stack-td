#!/bin/bash
# Refine a list of SAC files by a number of criteria for creating a list of
# stations when stacking

usage() {
	cat <<-END > $([ "$1" ] && echo "/dev/stdout" || echo "/dev/stderr")
	Usage: $(basename $0) (options) (filters) [list of SAC files]
	Options:
	   -h, --help           : Show help
	   -q                   : Do not quote the filenames.  This is done by
	                          default in case there are path separators,
	                          which are read as formatting by Fortran programs.
	Filters:
	   -d [delta1] [delta2] : Distance range
	   -a [az1] [az2]       : Azimuth range
	   -b [baz1] baz2]      : Backazimuth range
	   -lon [lon1] [lon2]   : Longitude
	   -lat [lat1] [lat2]   : Latitude
	END
	exit $([ "$1" ] || echo 1)
}

# Defaults
d1=0;      d2=360
az1=-360;  az2=720
baz1=-360; baz2=720
lon1=-360; lon2=360
lat1=-90;  lat2=90

[ $# -eq 0 ] && usage
while [ "$1" ]; do
	case "$1" in
		# Options
		-h|--help) usage 1;;
		-q) noquote=1; shift;;
		# Filters
		-d) d1="$2"; d2="$3"; shift 3;;
		-a) az1="$2"; az2="$3"; shift 3;;
		-b) baz1="$2"; baz2="$3"; shift 3;;
		-lon) lon1="$2"; lon2="$3"; shift 3;;
		-lat) lat1="$2"; lat2="$3"; shift 3;;
		*) [ -f "$1" ] && break || usage
	esac
done

saclst gcarc az baz stlo stla f "$@" |
awk -v noquote="$noquote" \
	-v d1="$d1"     -v d2="$d2" \
	-v az1="$az1"   -v az2="$az2" \
	-v baz1="$baz1" -v baz2="$baz2" \
	-v lon1="$lon1" -v lon2="$lon2" \
	-v lat1="$lat1" -v lat2="$lat2" '
	$2 >= d1 && $2 <= d2 && \
	$3 >= az1 && $3 <= az2 && \
	$4 >= baz1 && $4 <= baz2 && \
	$5 >= lon1 && $5 <= lon2 && \
	$6 >= lat1 && $6 <= lat2 {
		if (!noquote) $1 = "\""$1"\""
		print $1
	}'

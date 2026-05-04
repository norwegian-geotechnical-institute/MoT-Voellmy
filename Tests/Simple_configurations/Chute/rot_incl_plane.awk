# Create 10×20 grid with inclined plane rotated away from coordinate directions:

BEGIN {
	pi = 3.14159
	x0 = 0.0
	y0 = 0.0
	dx = 5.0
	nod = -9999.0
	nr = 10
	nc = 20
	mx = sin(40.0/180.0*pi)/cos(40.0/180.0*pi)
	my = sin(pi/6.0)/cos(pi/6.0)
	print("ncols       ", nc)
	print("nrows       ", nr)
	printf("xllcorner    %.2f\n", x0)
	printf("yllcorner    %.2f\n", y0)
	printf("cellsize     %.2f\n", dx)
	print("NODATA_value", nod)

	for (i=0; i<nr; i++) {
		for (j=0; j<nc-1; j++)
			printf("%6.2f  ", ((i+0.5)*my + (nc-j-0.5)*mx) * dx)
		printf("%6.2f\n", ((i+0.5)*my + 0.5*mx) * dx)
	}
}
# ' > rot_incl_plane.asc


################################################################

# Create initial condition:

{   if (NR < 7) print ;
    if (NR > 6) {
        for (i=1; i<=NF; i++) $i = ($i > 100.0 ? 1.0 ; 0.0) ;
        print
}
# ' rot_incl_plane.asc > rot_incl_plane_h0.asc


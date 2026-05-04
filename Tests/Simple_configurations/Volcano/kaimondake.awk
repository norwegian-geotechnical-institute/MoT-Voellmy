# Create "volcano" with tan θ = 0.6, height 900 m on a 5 m grid:

gawk '
BEGIN {
	pi = 3.14159 ; x0 = -1001.25 ; y0 = -1001.25 ; dx = 5.0
	nod = -9999.0 ; nr = 1001 ; nc = 1001 ; m = 0.6 ; zm = 900.0

	print("ncols       ", nc)
	print("nrows       ", nr)
	printf("xllcorner    %.2f\n", x0)
	printf("yllcorner    %.2f\n", y0)
	printf("cellsize     %.2f\n", dx)
	print("NODATA_value", nod)

	for (i=1; i<=nr; i++) {
		for (j=1; j<nc; j++) {
			z = 900.0 - m * dx * sqrt((i-501)**2 + (j-501)**2)
			printf("%6.2f  ", z > 0.0 ? z : 0.0)
		}
		printf("  0.00\n")
	}
}
' > kaimondake_dtm.asc


################################################################

# Create initial condition – uniform snow depth in a radius of 200 m around summit:

gawk '
BEGIN {
    h0 = 1.5 ; c = 501 ; r = 40
}

{   if (NR < 7) print ;
    if (NR > 6) {
        for (i=1; i<NF; i++)
            printf("%3.1f ", (NR-6-c)**2 + (i-c)**2 < r**2 ? h0 : 0.0 )
        printf("0.0\n")
    }
}
' kaimondake_dtm.asc > kaimondake_h0.asc


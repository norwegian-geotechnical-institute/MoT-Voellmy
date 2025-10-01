BEGIN {
    dx = 1.0 ;
    L = 100.0*(1.0-0.5*sqrt(2.0)) ;
    H = L + 1000.0 ;
    B = 50.0*sqrt(2.0) 
}

{
    if (NR < 7) print
    if (NR > 6) {
        for (i=0; i<=1000; i++) printf("%6.1f  ", H-i)
        for (i=1; i<=B; i++) {
            printf("%6.1f  ", 100.0 - sqrt(1.0e4 - (B-i*dx)**2))
        }
        for (i=1; i < 1000; i++) printf("   0.0  ")
        printf("   0.0\n") 
    } 
}

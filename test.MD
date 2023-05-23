COTAN TESTS

Genes = 5000 Cells = 500 calculateCoex() parallel: Time difference of 3.333824 secs Vanilla: Time difference of 4.235882 secs cellsUniformClustering() Parallel: Time difference of 2.820773 mins Vanilla: Time difference of 2.782857 mins

# tabula-muris-heart \# cells 13646 \# genes 617

### CalculateCoex()

-   Parallel: Time difference of 6.824892 secs

-   Vanilla: Time difference of 9.794934 secs

    Parallel = 30% faster than Vanilla

### cellsUniformClustering(GDI = 1.4)

-   Parallel: Time difference of 16.0639 mins --\> \~16 GB RAM (some peaks to \~19 GB) \*\*

-   Vanilla: Time difference of 15.53055 mins --\> \~15 GB RAM (some peaks to \~17/18 GB) \*\*

    Parallel 3.5% slower than Vanilla (with \~1 GB more RAM)

### mergedUniformCellsClusters(GDI = 2.1)

-   Parallel: Time difference of 10.63535 mins --\> on average 16 range(14, 18) GB RAM (some peaks to \~19 GB) \*\*

-   Vanilla: Time difference of 11.65808 mins --\> on average 15 range(15, 18) GB RAM (some peaks to \~18 GB) \*\*

    Parallel 9% faster than Vanilla (with \~1 GB more RAM)

# peripheal-blood \# cells 18306 \# genes 3658

### CalculateCoex()

-   Parallel: Time difference of 2.936496 mins

-   Vanilla: Time difference of 5.444983 mins

    Parallel 46% faster than Vanilla

### cellsUniformClustering(GDI = 2.0)

-   Parallel: Time difference of 11.74619 mins --\> \~16/16.5 GB RAM (some peaks to \~19 GB) \*\*

-   Vanilla: Time difference of 18.66191 mins --\> \~14.5 GB RAM (some peaks to \~17 GB) \*\*

    Parallel 37%. faster than Vanilla (with \~2/2.5 GB more RAM)

### mergedUniformCellsClusters(GDI=1.4)

-   Parallel: Time difference of 6.124006 mins --\> ~14/15 GB RAM (some peaks to \~18 GB) \*\*
-   Vanilla: Time difference of 13.10064 mins --\> ~13/14 GB RAM (some peaks to \~16 GB) \*\*
    
    Parallel 53% faster than Vanilla (with \~1 GB more RAM)

\*\* \[Compresa la memoria usata dal sistema\]
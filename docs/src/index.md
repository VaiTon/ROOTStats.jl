# Features

## ``\chi^2`` test

```@docs
chi2test
```

```@example
using FHist, ROOTStats

d1 = rand(10)
d2 = rand(10)

h1 = Hist1D(d1; binedges=0:0.01:1)
h2 = Hist1D(d2; binedges=0:0.01:1)

χ², ndf = chi2test(h1, h2, :WW)

print("chi2 = $χ², ndf = $ndf")
```


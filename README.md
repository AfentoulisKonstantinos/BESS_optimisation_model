# BESS Optimisation Model

The optimisation model maximises the profits of a BESS providing arbitrage, while participating in three different energy markets. 
The degradation cost of the BESS is 100£ per cycle. Thus, the degradation cost is not considered in the objective function since if considered, the annual cycles of the BESS are significantly lower, driving to degradation due to the maximum lifetime of ten years. For that reason, the degradation of the BESS before the end of the lifetime is avoided with a constraint that ensures that for the three year optimisation horizon the total cycles are lower than 1500.

## Requirements

- Install [Julia v1.10.4](https://julialang.org/downloads/platform/)
- Install Julia packages:
    ```
    using Pkg
    Pkg.add("CSV")
    Pkg.add("DataFrames")
    Pkg.add("Dates")
    Pkg.add("JuMP")
    Pkg.add("HiGHS")
    ```

## Run

``` 
git clone https://github.com/AfentoulisKonstantinos/BESS_optimisation_model
cd BESS_optimisation_model
julia BESS_optimiisation.jl
```

## Contact
Afentoulis Konstantinos: kwstasafentoulis@hotmail.com

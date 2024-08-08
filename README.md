# Calculation of couplexes from PICO experiments

- compatible with QIAcuity Software Suite 2.5.0.1 only
- sample names containing NTC will be removed
- supported Nanoplate formats 8.5k (use 13µl master mix) and 26k (use 42µl master mix)
- clusters with 0 counts are removed


*more information on usage will follow*

[App on shinyapps.io](https://thundert.shinyapps.io/calculate_couplexes/)


---

### How the PICO class is initialized, what its functions do and how the data frames are handled
```mermaid
graph TD;
%% classDef sub opacity:0.5
%% classDef note stroke:none
    subgraph init ["__init__()"]
        A(Upload and read csv)-->|self._calculate_clusters|B(self.df_clusters);
        B-->|self._general_formatting|C(self.df_clusters_formatted);
        C-->|self._general_filtering|D(self.df_filtered_prelim);
        D-->|self._format_for_lambda_hist|E(self.df_lambda);
        D-->|self._calculate_couplexes|G(self.df_couplexes);
    end
        E-->|self.get_lambda_range|F(histogram of overall\nλ range for sidebar);
        G-->|self.get_couplex_plot|H(violin plot of couplexes\nfor main panel)
        J-->|self.get_couplex_plot|H
        H-->|download|K(.pdf)
        G-->|download|L(.csv)
        J-->|download|M(.csv)
    subgraph lambda ["if input.lambda_filter()"]
        G-->|self.lambda_filtering|J(self.df_couplexes_filtered)
        I(filter values from silder)-->J
    end
    G-->|self.get_lambda_ranges|N(λ range plot for\nmain panel)
    J-->|self.get_lambda_ranges|N
    N-->|download|O(.pdf)


    %% class noteA note   
```

### to dos
- ```cluster_calculation.py``` still runs with pandas, while the rest runs with polars. This necessiates the conversion of dataframe types at some points.
- 
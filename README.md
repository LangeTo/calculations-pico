# Calculation of couplexes from PICO experiments

compatible with QIAcuity Software Suite 2.5.0.1

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
        D-->|self._format_for_lambda|E(self.df_lambda);
        D-->|self._calculate_couplexes|G(self.df_couplexes);
    end
        E-->|self.get_lambda_range|F(histogram of overall\nλ range for sidebar);
        G-->|self.get_couplex_plot|H(violin plot of couplexes\n for main panel)
        J-->|self.get_couplex_plot|H
        H-->|download|K(.pdf)
        G-->|download|L(.csv)
        J-->|download|M(.csv)
    subgraph lambda ["if input.lambda_filter()"]
        G-->|self.lambda_filtering|J(self.df_couplexes_filtered)
        I(filter values from silder)-->J
    end

    %% subgraph subinit [" "]
    %%     init
    %%     noteA([I AM THE FIRST NOTE])
    %% end

    %% class noteA note
    
```

<!-- (re)deploying the app

```powershell
rsconnect deploy shiny C:\Users\tl100\PycharmProjects\shiny_amulator --name thundert --title calculate_couplexes
``` -->

<!-- rename git repository
* disconnect from remote ```git remote rm origin```
* add the new remote branch ```git remote add origin https://github.com/LangeTo/calculations-pico.git```
* then set upstream branch ```git push --set-upstream origin master``` -->

<!-- might also be interesting instead of shiny: https://docs.bokeh.org/en/latest/index.html#  -->
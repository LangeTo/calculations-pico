# Calculation of couplexes from PICO experiments

compatible with QIAcuity Software Suite 2.5.0.1

*more information on usage will follow*

[App on shinyapps.io](https://thundert.shinyapps.io/calculate_couplexes/)


---
Issues with deploying the app to 

downgrade to python 3.11.9 necessary
```powershell
rsconnect write-manifest shiny . 
```
Then, still remove version of ```mkl-service``` and comment ```pywin32```. However, I think these packages are not used anyhow

```powershell
pip freeze > requirements.txt
```

```powershell
conda list -e > requirements.txt
```

(re)deploying the app

```powershell
rsconnect deploy shiny C:\Users\tl100\PycharmProjects\shiny_amulator --name thundert --title calculate_couplexes
```


<!-- rename git repository
* disconnect from remote ```git remote rm origin```
* add the new remote branch ```git remote add origin https://github.com/LangeTo/calculations-pico.git```
* then set upstream branch ```git push --set-upstream origin master``` -->
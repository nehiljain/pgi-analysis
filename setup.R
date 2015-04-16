list_of_packages <- c("formatR", "highr", "markdown", "knitr","crayon", "jsonlite", "rstudioapi","httr","dplyr", "DBI","mailR","xlsx")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')

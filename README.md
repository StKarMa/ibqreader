## hello person!

thanks for being here.  
this might be an accident or you are one of the few persons in the world  
using isobarquant for protein quantification.

this little package will help to glimpse your data and maybe even do other stuff
Its written in R and therefore you should have R (studio) running.

copy/paste the the next lines into R to install the development version 

> library(devtools)  
> install_github("StKarMa/ibqreader@mshguest")  
> library(ibqreader)

to get a first idea what we can do, we can put some example data into our working directory 


> path2data <- paste0(system.file(package = "ibqreader"), "/extdata/example_data")
> dir.create("example_data")
> file.copy(path2data, "example_data", recursive=TRUE)








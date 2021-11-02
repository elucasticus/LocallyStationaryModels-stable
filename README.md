# **LocallyStationaryModels**
***authors:*** Luca Crippa, Giacomo De Carlo<br>
***mailto:*** <luca7.crippa@mail.polimi.it>, <giacomo.decarlo@mail.polimi.it><br>
## **Latest updates**
*01/11* Added the function `cv.lsm(...)` to calculate the mean squared error via cross-validation.
*30/10* ⚠️ The use of the function `variogramlsm(...)` is now deprecated. Use `variogram.lsm(...)` instead.
## **Installation**
### **Step 1: install R and Rstudio (optional)**
#### **Windows:**
Download R from <https://cran.r-project.org/bin/windows/base/> and (optional) Rstudio from <https://www.rstudio.com/products/rstudio/download/>.
#### **Arch Linux:**
Open the terminal then type

    sudo pacman -S r

and install the `r` package for R.
To get Rstudio (optional), install an AUR helper of your choice such as <https://github.com/Jguer/yay> and then on the terminal

    yay -S rstudio-desktop-bin
Refer to <https://wiki.archlinux.org/title/r> for further details.
#### **Ubuntu:**
Refer to <https://cran.r-project.org/bin/linux/ubuntu/> for R.
For Rstudio (optional) you can download the `.deb` file from <https://www.rstudio.com/products/rstudio/download/>.
### **Step 2: install devtools and LocallyStationaryModels**
Open R then type

    install.packages("devtools")
    library(devtools)
    devtools::install_github("giacomodecarlo/LocallyStationaryModels")

`LocallyStaionaryModels` requires `ggplot2`, `ggforce` and `cowplot` to plot the results. If needed you can install them from CRAN via the commands

    install.packages("ggplot2")
    install.packages("ggforce")
    install.packages("cowplot")

You can access the dataset *meuse* to run the examples via the `sp` library. Just open R and type

    install.packages("sp")
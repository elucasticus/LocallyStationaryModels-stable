# **LocallyStationaryModels**
***authors:*** Luca Crippa, Giacomo De Carlo<br>
***mailto:*** <luca7.crippa@mail.polimi.it>, <giacomo.decarlo@mail.polimi.it><br>
## **Installation**
### **Step 1: install R and Rstudio (optional)**
#### **Windows:**
Download R from <https://cran.r-project.org/bin/windows/base/> and (optional) Rstudio from <https://www.rstudio.com/products/rstudio/download/>.
#### **Arch Linux:**
Open the terminal then type

    sudo pacman -S r

to install the `r` package.
To install Rstudio (optional) install an AUR helper of your choice such as yay <https://github.com/Jguer/yay> then on the terminal

    yay -S rstudio-desktop-bin
#### **Ubuntu:**
Refer to <https://cran.r-project.org/bin/linux/ubuntu/> for `r`.<br>
For Rstudio (optional) you can download the `.deb` file from <https://www.rstudio.com/products/rstudio/download/>.
### **Step 2: install devtools and LocallyStationaryModels**
Open `r` then type

    install.packages("devtools")
    library(devtools)
    devtools::install_github("giacomodecarlo/LocallyStationaryModels")

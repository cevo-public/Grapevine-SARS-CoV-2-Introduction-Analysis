FROM rocker/tidyverse:4.0.3
WORKDIR /app/

# Install IQTree
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.6/iqtree-2.0.6-Linux.tar.gz
RUN tar -xvf iqtree-2.0.6-Linux.tar.gz

# Install R packages
COPY install_packages.R .
RUN Rscript /app/install_packages.R

# Install Python
RUN apt-get update && apt-get install -y python3-pip
COPY database/python/requirements.txt database/python/requirements.txt
RUN pip3 install -r database/python/requirements.txt

COPY . .
RUN chmod +x /app/main.sh

# Start program
ENTRYPOINT cd /app && bash /app/main.sh

# Contributing to BioTEA-box
The box is made from a base image that contains all required packages, the `biotea-base` image. The code is inserted from the base image in the container to allow the separate management of dependencies and code. In this way, the required R libraries are frozen, and automated tests (such as here on github) do not take three hours of compile time to complete. 

The box code is divided in modules. Each module is called in a particular way by the entrypoint (`entrypoint.R`), so please read the entrypoint file first. A template for new modules is inside `src/box/modules/template/`, so you can start wih that if you want to make a new module.

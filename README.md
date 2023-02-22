![](draft/assets/images/repo-banner.png)

# Kerswell et al. (2023; G3)

This repository provides all materials for the manuscript *Computing Rates and Distributions of Rock Recovery in Subduction Zones* (Kerswell et al., 2023; G3). You can find the preprint [here](https://buchanankerswell.com/assets/pdf/kerswell-et-al-marx-g3-2023.pdf).

This repository includes:

- An R script to download all data required to compile the study
- R scripts to run the analyses and reproduce all results and figures
- A Makefile and run.sh script to easily compile the study
- The complete manuscript written in Rmarkdown

This repository is self-contained but requires the following software (all open-source).

## Prerequisite software

### R

This study is written in R. Follow the instructions at [R's homepage](https://www.r-project.org) to download and install the latest release of R on your machine.

## Running the study

```
# Clone this repository
git clone https://github.com/buchanankerswell/kerswell_et_al_marx.git

# Change into the directory
cd kerswell_et_al_marx

# Use Makefile to compile
make
```

This will check for required R packages and try to install missing packages automatically.

If all packages are found and available it will proceed to run the study with some initial prompts from the user. The study takes about 2.5-3hr to run on my MacBook Air (M1 8GB, 2020) with 8 cores computing in parallel and 20 jackknife samples. The study takes about 7 minutes to run without resampling.

# Open Science Framework

This repository can also be found at the official [OSF repository](https://osf.io/3emwf/).

# Funding

This project was supported by the NSF grant OISE 1545903 to M. Kohn, S. Penniston-Dorland, and M. Feineman

# License

MIT License

Copyright (c) 2021 Buchanan Kerswell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


HTML = index.html
# Get name of all Rmd files
RMD = index.Rmd
# local assets
lASSETS := $(shell find assets)
# global assets
gASSETS := $(shell find ../assets)

# assets
$(HTML): $(RMD) $(lASSETS) $(gASSETS)
	Rscript -e "rmarkdown::render('$(RMD)', output_file = '$(HTML)')"
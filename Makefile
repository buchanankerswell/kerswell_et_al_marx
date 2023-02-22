R = R/classify.R R/download-data.R R/functions.R R/marx-vis.R R/packages.R R/preprocess.R R/summary.R R/trace-grids.R R/trace-marx.R
DATAPURGE = data/classification-stats-summary-*.csv data/preprocessed.RData data/marx_classified_*/*.csv data/marx-class-summary-*.RData data.zip
DATACLEAN = data draft/assets/r/*.csv draft/assets/r/*.RData draft/assets/r/marx_classified_*
FIGSPURGE = figs draft/assets/figs

all: $(R)
	@./run.sh

purge:
	@rm -rf $(DATAPURGE) $(FIGSPURGE)

clean: purge
	@rm -rf $(DATACLEAN)

.PHONY: all purge clean
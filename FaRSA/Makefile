# SUBDIRS = src tests problems exes
SUBDIRS = src tests functions ftests

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for d in $(SUBDIRS); \
	do \
	$(MAKE) --directory=$$d clean; \
	done

veryclean:
	for d in $(SUBDIRS); \
	do \
	$(MAKE) --directory=$$d veryclean; \
	done

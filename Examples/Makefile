PROJECT_DIRECTORIES := $(sort $(dir $(wildcard */)))
ZIP_NAMES := $(patsubst %/, %.zip, $(PROJECT_DIRECTORIES))
all:zip
test:
	echo $(ZIP_NAMES)
$(ZIP_NAMES):%.zip:%
	@rm -f $@ 
	@zip -r $@ $<
zip:$(ZIP_NAMES)

clean:
	rm -f  $(ZIP_NAMES)

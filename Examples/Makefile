PROJECT_DIRECTORIES := $(sort $(dir $(wildcard */)))
ZIP_NAMES := $(patsubst %/, %.zip, $(PROJECT_DIRECTORIES))
all:zip
test:
	echo $(ZIP_NAMES)
$(ZIP_NAMES):%.zip:%
	@rm -f $@ 
	(cd $< && zip -r ../$@ *)
all_zip:
	zip -r Examples.zip */
zip:$(ZIP_NAMES) all_zip

clean:
	rm -f  $(ZIP_NAMES)

PROJECT_DIRECTORIES := $(sort $(dir $(wildcard */)))
ZIP_NAMES := $(patsubst %/, %.zip, $(PROJECT_DIRECTORIES))
all:zip
test:
	echo $(ZIP_NAMES)
$(ZIP_NAMES):%.zip:%
	zip -r $@ $<
zip:$(ZIP_NAMES)

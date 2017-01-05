#######################
# DEVELOPER UTILITIES #
#######################

new-clone-init:
	git submodule update --init --recursive

update-adapter-clipper:
	cd adapter-clipper && git pull origin develop && cd ..

update-dge-analysis:
	cd dge-analysis && git pull origin master && cd ..

update-read_counter:
	cd read_counter && git pull origin master && cd ..

##################
# USER UTILITIES #
##################

clean-test:
	find test/f* -maxdepth 1 -mindepth 1 -type d ! -name 'raw' | xargs rm -rf
	rm -rf test/repro-archive/

test-run:
	cd test && sh pipeline_testrun.sh && cd ..

unbuild-test-reference:
	find test/ref* -maxdepth 1 -mindepth 1 ! -name '*.fa.gz' | xargs rm -rf

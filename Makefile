#######################
# DEVELOPER UTILITIES #
#######################

new-clone-init:
	git submodule update --init --recursive

update-dge-analysis:
	cd dge-analysis && git pull origin master && cd ..

update-read_counter:
	cd read_counter && git pull origin master && cd ..

##################
# USER UTILITIES #
##################

clean-test:
	rm -rf $(find test/f* -depth 1 -type d ! -name 'raw')

test-run:
	cd test && sh pipeline_testrun.sh && cd ..

unbuild-test-reference:
	rm -rf $(find test/reference -depth 1 ! -name '*.fa.gz')
	

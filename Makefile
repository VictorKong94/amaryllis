new-clone-init:
	git submodule update --init --recursive

update-dge-analysis:
	cd dge-analysis && git pull origin master && cd ..

update-read_counter:
	cd read_counter && git pull origin master && cd ..

test-run:
	cd test && sh pipeline_testrun.sh && cd ..

new:
	git submodule update --init --recursive

update-clipper:
	cd adapter-clipper && git pull origin develop && cd ..

update-dge:
	cd dge-analysis && git pull origin master && cd ..

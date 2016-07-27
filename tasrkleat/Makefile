build:
	# build the docker image
	@docker build -t ${USER}/tasrkleat . 2>&1 | tee build.log

run: build
	# Run docker container for debugging purpose locally
	@docker run --rm -ti \
	    -v ${PWD}/experiment:/experiment \
	    -v ${PWD}/app:/app \
	    ${USER}/tasrkleat  /bin/bash

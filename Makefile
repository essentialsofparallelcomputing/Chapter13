All: ShallowWater

ShallowWater:
	cd OpenACC/ShallowWater/ && mkdir build && cd build && export GRAPHICS_TYPE=JPEG && \
	   cmake -DENABLE_GRAPHICS=1 .. && make && ../run.sh

clean:
	cd OpenACC/ShallowWater && rm -rf build

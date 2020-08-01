All: ShallowWater

ShallowWater:
	cd OpenACC/ShallowWater/ && mkdir build && cd build && cmake .. && make && ../run.sh

clean:
	cd OpenACC/ShallowWater && rm -rf build

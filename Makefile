build:
	docker build -t orchestra .

build-singularity:
	singularity build --fakeroot orchestra.sif Singularity.def

run-singularity:
	singularity run orchestra.sif

shell-singularity:
	singularity shell orchestra.sif

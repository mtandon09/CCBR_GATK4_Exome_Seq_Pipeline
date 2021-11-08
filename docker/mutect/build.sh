# Build image
docker build -f Dockerfile --tag=ccbr_mutect:v0.1.0 .

# Tag image with version and reset latest
docker tag ccbr_mutect:v0.1.0 skchronicles/ccbr_mutect:v0.1.0
docker tag ccbr_mutect:v0.1.0 skchronicles/ccbr_mutect
docker tag ccbr_mutect:v0.1.0 nciccbr/ccbr_mutect:v0.1.0
docker tag ccbr_mutect:v0.1.0 nciccbr/ccbr_mutect

# Push image to DockerHub
docker push skchronicles/ccbr_mutect:v0.1.0
docker push skchronicles/ccbr_mutect:latest
docker push nciccbr/ccbr_mutect:v0.1.0
docker push nciccbr/ccbr_mutect:latest

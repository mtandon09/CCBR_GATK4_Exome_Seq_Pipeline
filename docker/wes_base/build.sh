# Build image
docker build -f Dockerfile --tag=ccbr_wes_base:v0.1.0 .

# Tag image with version and reset latest
docker tag ccbr_wes_base:v0.1.0 skchronicles/ccbr_wes_base:v0.1.0
docker tag ccbr_wes_base:v0.1.0 skchronicles/ccbr_wes_base
docker tag ccbr_wes_base:v0.1.0 nciccbr/ccbr_wes_base:v0.1.0
docker tag ccbr_wes_base:v0.1.0 nciccbr/ccbr_wes_base

# Push image to DockerHub
docker push skchronicles/ccbr_wes_base:v0.1.0
docker push skchronicles/ccbr_wes_base:latest
docker push nciccbr/ccbr_wes_base:v0.1.0
docker push nciccbr/ccbr_wes_base:latest

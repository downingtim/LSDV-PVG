There are two docker files for building pggb and dependancies. The dockerfile Dockerfile.pggb is based on a pre-built image that worked well on our test data, whereas some later pggb version failed to properly segment some test data. If the latest version of pggb needs to be incorporated, Dockerfile.pggb.latest can be used. 

To build the image cleanly, you can use the command
docker build --no-cache -t chandanatpi/tpi:pggbV01 -f docker/Dockerfile.pggb.latest . 
from the base directory of the package.

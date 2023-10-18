1. Update `vcfdist` version in `globals.h`
2. Replace `vcfdist --version` text in `src/README.md`
3. Update `vcfdist` version in `README.md`
4. Build and deploy new Docker image
```
cd ~/vcfdist/src
sudo docker login -u timd1
sudo docker build --no-cache -t vcfdist .
sudo docker image tag timd1/vcfdist:latest timd1/vcfdist/vX.X.X
sudo docker image push timd1/vcfdist:vX.X.X
```
5. Make final Git commit of changes
6. Make new release and tag on Github

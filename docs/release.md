1. Update `vcfdist` version in `globals.h`
2. Replace `vcfdist --version` text in `src/README.md`
3. Update `vcfdist` version in `README.md`
4. Make final Git commit of changes
```bash
git commit -m "version bump (vX.X.X)"
git checkout -b master
git merge dev
git push origin master
git checkout -b dev
```
5. Make new release and tag on Github
6. Build and deploy new Docker image
```bash
cd ~/vcfdist/src
sudo docker login -u timd1
sudo docker build --no-cache -t vcfdist .
sudo docker image tag timd1/vcfdist:latest timd1/vcfdist/vX.X.X
sudo docker image push timd1/vcfdist:vX.X.X
```

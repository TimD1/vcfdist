## Release Notes
More information on previous releases of vcfdist can be found [here](https://github.com/timd1/vcfdist/releases).

## Steps to Create New Release

1. Compare against previous releases using the `run` and `demo` scripts, updating `output.txt` for `demo` if necessary

2. Update `vcfdist` version in `src/globals.h`

3. Update `vcfdist` version in `README.md`

4. Cache vcfdist wiki.
```
cd vcfdist/docs
git clone https://github.com/TimD1/vcfdist.wiki.git vX.X.X
rm -rf vX.X.X/.git
```

5. Make final Git commit of changes
```bash
git add src/globals.h README.md docs/
git commit -m "version bump (vX.X.X)"
git checkout master
git merge dev
git push origin master
git checkout dev
```

6. Make new release and tag on Github

7. Build and deploy new Docker image
```bash
cd ~/vcfdist/src
sudo docker login -u timd1
sudo docker build --no-cache -t timd1/vcfdist .
sudo docker image ls
sudo docker image tag timd1/vcfdist:latest timd1/vcfdist:vX.X.X
sudo docker image push timd1/vcfdist:vX.X.X
sudo docker image push timd1/vcfdist:latest
```
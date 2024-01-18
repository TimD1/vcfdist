1. Benchmark against `run` and `demo` scripts
- update `output.txt` for `demo` if necessary

2. Update `vcfdist` version in `globals.h`

3. Update `vcfdist` help text in `src/README.md`
```bash
cd src
make
./vcfdist >> README.md
```
Manually delete previous text

4. Update `vcfdist` version in `README.md`

5. Make final Git commit of changes
```bash
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

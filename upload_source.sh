rm -rf ./source
mkdir source
find ./src/ -name "*.pyc" -print0 | xargs -0 rm -rf
cp -R ./src ./source/
tar czf code.tar.gz source
scp code.tar.gz hto@aludra.usc.edu:/home/scf-11/hto/geocrowd

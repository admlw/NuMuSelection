for f in `cat test.list`; do 
  root -l -q ${f} test.C;
done;

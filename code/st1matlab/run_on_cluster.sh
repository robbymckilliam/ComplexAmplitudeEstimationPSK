#!/bin/bash
filename=$1
for file in $filename 
do
  if [ -f "$file" ]
  then
    tmp=${file/run_}
    wrapper=${tmp/.sh/_wrap.sh}
    matfile=${tmp/.sh} 
    if [ -f "$wrapper" ]
    then
      echo "Skipping $wrapper as file already exists."
    else
      echo "Processing $file"
      echo "    Generating $wrapper..."
      echo "#!/bin/bash" >> $wrapper
      echo ./$file /usr/local/MATLAB/MATLAB_Compiler_Runtime/v80 \"$matfile\" >> $wrapper
      chmod u+x $wrapper
      echo "    Executing qsub $wrapper"
      qrun `qsub -d /home/polloka/st1matlab/standalone/ $wrapper` 
    fi
  else
    echo "$file does not exist."
  fi
done
sleep 2s
echo
echo "QSTAT:"
echo
qstat

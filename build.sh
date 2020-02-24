#!/bin/csh

if ( !($?COATJAVA) ) then
  set my_java = "6.3.1"
  set java_path = "/group/clas12/packages/coatjava"
  setenv COATJAVA ${java_path}/${my_java}/
  echo "COATJAVA version is set to the defaul $my_java. Update this script, if want a newer version. "
else
  echo "COATJAVA is set to \$COATJAVA"
endif

setenv MALLOC_ARENA_MAX 1

if ( $#argv == 0 ) then
  set my_script = "ana_2p2"
else if ($#argv == 1) then
  if ($1 == "-help") then
   echo "To compile all java codes, execute "./build.sh" "
   echo "To compile name.java, execute "build.sh name" "
   exit 0
  else
   set my_script = "$1"
  endif
else
  echo "Too many arguments"
endif
  echo "Compiling java script $my_script"
  javac -cp "$COATJAVA/lib/clas/*:$COATJAVA/lib/utils/*:." ${my_script}.java

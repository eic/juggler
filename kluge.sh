#!/bin/bash
# - KLUGE: the `sed -i ... $klugedFile` commands work around
#   Gaudi `listcomponent`'s preference of "installed" libraries
#   (maybe `$LD_LIBRARY_PATH`) over those we have just built; the regexp
#   forces `listcomponents` to use the newly built library in the `pwd`; since
#   this will write unecessary `./` characters to the output `*.components` file,
#   a second regexp is used to remove them
for klugedFile in `find build -name "build.make" -print | grep -E 'Jug.*Plugins\.dir'`; do
  echo "KLUGE $klugedFile"
  #mv $klugedFile{.bak,} # restore backups (uncomment this if you are brave enough to edit this script)
  components=$(grep listcomponents $klugedFile | sed 's/^.*--output //g' | awk '{print $1}') # components file
  echo "COMPONENTS $components"
  sed -i.bak '/listcomponents/{s/ libJug/ .\/libJug/;}' $klugedFile # make definite relative path
  sed -i     's/listcomponents.*$/& \&\& sed -i "s\/\\\.\\\/\/\/g" '${components}'/g' $klugedFile # remove ./ from $components
  #echo "DIFF:"; diff $klugedFile{.bak,} # show what the kluge did
done

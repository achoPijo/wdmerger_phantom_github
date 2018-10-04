#!/bin/sh

# The version filename
VF=../bin/phantom_version

#Important details about the commit
TimeNow=$(date)
GitDate=$(git log -n 1 | grep Date | cut -c6-)
GitId=$(git log -n 1 --pretty=oneline | cut -c1-41)
IsModified=$(git diff-index --name-only HEAD --)
GitNow=$(git describe --tag)
test -z "$IsModified" || IsMod="true"

# Write to file
echo "Phantom version details:" > $VF
echo "Date this version was committed to the repository: $GitDate" >> $VF
echo "Date Phantom was compiled: $TimeNow" >> $VF
echo "Commit tag:     $GitNow" >> $VF
echo "Full commit ID: $GitId" >> $VF
if IsMod=="true"
   then
      echo "Locally modified files:" >> $VF;
      echo "$IsModified" >> $VF;
fi



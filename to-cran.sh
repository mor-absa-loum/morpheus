#!/bin/bash

# Synchronize package folder
rsync -a --delete pkg/ pkg-cran/
cd pkg-cran

# Replace all used greek letters by beta, mu, ...etc
for file in `find -name *.R`; do
	if [ -f $file ]; then
		sed -i 's/μ/mu/g' $file
		sed -i 's/β/beta/g' $file
		sed -i 's/λ/lambda/g' $file
		sed -i 's/Σ/Sigma/g' $file
	fi
done

if [ "$1" == "i" ]; then
	# Install package
	echo "rp()" | R --slave
fi

cd ..

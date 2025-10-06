#!/bin/bash

# Set the number of subjects and surrogate runs
SUBS=30
NSURR=1000

read -p "Do you want to download the required MATLAB packages? (Brain Connectivity Toolbox, parc_plotter) (y/n) " download_packages

if [ "$download_packages" == "y" ]; then
    # Create the target directory if it doesn't exist
    mkdir -p funcs
    # Download package from GitHub
    git clone https://github.com/faskowit/parc_plotter.git funcs/parc_plotter

    # Download BCT
    FILE_ID="1DmMvRnferBfGe057O-sZwB5jL4j8w1Hu"
    FILE_NAME="BCT.zip"
    wget --no-check-certificate 'https://docs.google.com/uc?export=download&id='"$FILE_ID" -O "$FILE_NAME"
    unzip funcs/${FILE_NAME} -d funcs
    rm funcs/${FILE_NAME}
else
    echo "Warning: If the required functions are not available and appropriate changes are not made in the code to handle it, the functions will throw an error."
fi

# Loop through each subject
for SUB in $(seq 1 $SUBS); do
    # Run the code to obtain empirical (observed) results
    matlab -nodisplay -nosplash -r "sub=$SUB; observed/observed"

    # Run the surrogate code for each surrogate run
    for SURR in $(seq 1 $NSURR); do
        matlab -nodisplay -nosplash -r "sub=$SUB; surr=$SURR; surrogate/surrogate"
    done

    # Run the surrogate correction code
    matlab -nodisplay -nosplash -r "sub=$SUB; nsurr=$NSURR; surr_correction"
done

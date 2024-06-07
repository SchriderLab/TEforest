# this script splits training data based on chromosome.
 
while getopts s: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
    esac
done

echo "Current directory: $(pwd)"
featvec_dir="$(pwd)/feature_vectors_het_reference"
training_folder=$(pwd)/2L2R_reference
validation_folder=$(pwd)/3L3RX_reference

[ -d $training_folder ] && echo "Output directory exists" || mkdir $training_folder
[ -d $validation_folder ] && echo "Output directory exists" || mkdir $validation_folder

for dir in "$featvec_dir/*" #[a-z][0-9][a-z][0-9]
    do echo $dir
    # Remove files with size 800 bytes. This means they have no data.
    find $dir -type f -size 800c -exec rm -f {} +
    find $dir -type f -name "*2L*" -exec cp {} $training_folder/ \; 
    find $dir -type f -name "*2R*" -exec cp {} $training_folder/ \; 
    touch $training_folder/$sample.txt
    find $dir -type f -name "*3L*" -exec cp {} $validation_folder/ \; 
    find $dir -type f -name "*3R*" -exec cp {} $validation_folder/ \; 
    find $dir -type f -name "*X*" -exec cp {} $validation_folder/ \; 
    touch $validation_folder/$sample.txt
done
